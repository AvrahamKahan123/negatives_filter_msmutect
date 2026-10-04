"""Loads config.json and makes the rest of the repo importable.

Every cluster path lives in config.json, so the scripts and submit files carry none. Any
setting can be overridden at submit time by an environment variable of the same name
upper-cased (condor's `getenv = true` passes your shell through to the jobs):

    WORK_ROOT=/storage/bfe_maruvka/avrahamk/scratch_run condor_submit step1_filter.sub
"""

import errno
import json
import os
import sys

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
CONFIG_FILE = os.path.join(THIS_DIR, "config.json")

with open(CONFIG_FILE) as _f:
    _settings = {k: v for k, v in json.load(_f).items() if not k.startswith("_")}


def setting(name, default=None):
    """config.json value, overridden by the upper-cased environment variable if set."""
    from_env = os.environ.get(name.upper())
    if from_env is not None:
        if isinstance(_settings.get(name), bool):
            return from_env.strip().lower() in ("1", "true", "yes", "on")
        return from_env
    return _settings.get(name, default)


INPUT_DIR = setting("input_dir")
INPUT_SUFFIX = setting("input_suffix", ".full.mut.tsv.gz")
WORK_ROOT = setting("work_root")
MSMUTECT_PATH = setting("msmutect_path")
GZIP_INTERMEDIATES = setting("gzip_intermediates", True)
DELETE_FILTERED_AFTER_GERMLINE = setting("delete_filtered_after_germline", False)

# repo_root defaults to the parent of this directory, which is correct whenever the jobs
# run from the checkout. Put it on sys.path so steps 2 and 3 can import the shared packages.
REPO_ROOT = setting("repo_root") or os.path.dirname(THIS_DIR)
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

# copies recovered from the backup host for samples the main storage could not read
# (fetch_one_unread.sh). They are byte-identical in shape to the originals in INPUT_DIR,
# so they feed step 1 exactly the same way -- see filter_unread.py.
UNREAD_DIR = os.path.join(WORK_ROOT, "unread")

# the three output directories named in Plan.md
FILTERED_DIR = os.path.join(WORK_ROOT, "filtered")
FILTERED_WGERMLINE_DIR = os.path.join(WORK_ROOT, "filtered_wgermline")
MSMUTECT_RESULTS_DIR = os.path.join(WORK_ROOT, "msmutect_results")

JOBS_DIR = os.path.join(THIS_DIR, "jobs")
LOGS_DIR = os.path.join(WORK_ROOT, "logs")

# intermediates keep the input's shape; only the compression may differ
INTERMEDIATE_SUFFIX = ".full.mut.tsv.gz" if GZIP_INTERMEDIATES else ".full.mut.tsv"


def source_path(sample: str) -> str:
    return os.path.join(INPUT_DIR, sample + INPUT_SUFFIX)


def filtered_path(sample: str) -> str:
    return os.path.join(FILTERED_DIR, sample + INTERMEDIATE_SUFFIX)


def filtered_wgermline_path(sample: str) -> str:
    return os.path.join(FILTERED_WGERMLINE_DIR, sample + INTERMEDIATE_SUFFIX)


def msmutect_output_prefix(sample: str) -> str:
    return os.path.join(MSMUTECT_RESULTS_DIR, sample)


def msmutect_output_path(sample: str) -> str:
    return msmutect_output_prefix(sample) + ".full.mut.tsv"


def sample_name(path: str) -> str:
    """TCGA-W5-AA2K.full.mut.tsv[.gz] -> TCGA-W5-AA2K"""
    name = os.path.basename(path)
    if name.endswith(".gz"):
        name = name[: -len(".gz")]
    for suffix in (".filtered.full.mut.tsv", ".full.mut.tsv"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return os.path.splitext(name)[0]


OK, MISSING, UNREADABLE = "ok", "missing", "unreadable"


def path_state(path: str) -> tuple:
    """-> (OK | MISSING | UNREADABLE, detail or None).

    os.path.exists() returns False for BOTH an absent file and one the filesystem cannot
    stat -- /storage throwing EIO, a dead mount, a permissions problem. Those mean opposite
    things here: absent is normal (the previous step has not run yet), unreadable is a
    storage fault that would otherwise drop the sample from the batch without a word.
    """
    try:
        os.stat(path)
        return OK, None
    except FileNotFoundError:
        return MISSING, None
    except OSError as e:
        name = errno.errorcode.get(e.errno, str(e.errno))
        return UNREADABLE, f"{name}: {e.strerror}"


def require(path: str, what: str):
    if not path:
        raise SystemExit(f"ERROR: {what} is not set in {CONFIG_FILE}")
    if not os.path.exists(path):
        raise SystemExit(f"ERROR: {what} does not exist: {path}\n"
                         f"       Fix it in {CONFIG_FILE} (or override with the "
                         f"environment variable).")
    return path


def describe() -> str:
    return "\n".join([
        f"  input_dir            {INPUT_DIR}",
        f"  work_root            {WORK_ROOT}",
        f"  msmutect_path        {MSMUTECT_PATH}",
        f"  repo_root            {REPO_ROOT}",
        f"  gzip_intermediates   {GZIP_INTERMEDIATES}",
    ])


if __name__ == "__main__":
    print(f"config from {CONFIG_FILE}")
    print(describe())
