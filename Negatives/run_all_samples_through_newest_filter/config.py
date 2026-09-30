import os
import sys

# --- make results_postprocessing importable, the same way analyze_full_tcga_mut_set.py expects ---
# results_postprocessing/cancer_data.py does `from SamplesDB import ...` (bare, not package-relative),
# so results_postprocessing itself must be on sys.path in addition to the project root.
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(THIS_DIR, "..", ".."))
RESULTS_POSTPROCESSING_DIR = os.path.join(PROJECT_ROOT, "results_postprocessing")
for _p in (PROJECT_ROOT, RESULTS_POSTPROCESSING_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)

# --- samples under test ---
SAMPLES = [
    "TCGA-A6-5661",
    "TCGA-AJ-A3BH",
    "TCGA-AP-A05N",
    "TCGA-FI-A2D4",
    "TCGA-OR-A5LB",
]

# --- remote source (see Negatives/run_all_samples_through_newest_filter/Plan.md) ---
REMOTE_HOST = "tech-ui01.hep.technion.ac.il"
REMOTE_USER = "avrahamk"
REMOTE_RUN_DIR = (
    "/storage/bfe_maruvka/gaiafr/Research_project/WGS_SNVs_indels_analysis_project"
    "/MS_Analysis/Google_Cloud_MSMuTect_final_output_new_May2025_run"
)
REMOTE_FILT_DIR = f"{REMOTE_RUN_DIR}/MSMuTect_called_mut_filt_fixed"
REMOTE_FULL_DIR = REMOTE_RUN_DIR


def remote_filt_path(sample: str) -> str:
    return f"{REMOTE_FILT_DIR}/{sample}.called.filt.mut.tsv.gz"


def remote_full_path(sample: str) -> str:
    return f"{REMOTE_FULL_DIR}/{sample}.full.mut.tsv.gz"


# --- local layout (everything lives under this directory, per Plan.md) ---
DATA_DIR = os.path.join(THIS_DIR, "data")
FILT_DIR = os.path.join(DATA_DIR, "filt")          # downloaded *.called.filt.mut.tsv.gz
FULL_DIR = os.path.join(DATA_DIR, "full")          # downloaded *.full.mut.tsv.gz (one at a time)

RESULTS_DIR = os.path.join(THIS_DIR, "results")
OLD_YOSSI_DIR = os.path.join(RESULTS_DIR, "old_yossi_filtered")     # step 2 output
RERUN_DIR = os.path.join(RESULTS_DIR, "from_file_rerun")            # step 3 output (msmutect --from_file)
NEW_YOSSI_DIR = os.path.join(RESULTS_DIR, "new_yossi_filtered")     # yossi_filter applied to step 3 output
COMPARISON_DIR = os.path.join(RESULTS_DIR, "comparison")            # step 4 output

MSMUTECT_PATH = os.path.expanduser("~/MaruvkaLab/MSMuTect/msmutect.sh")

# locus columns carved out for --from_file (KNOWN_GERMLINE_VARIANT required by the normal
# histogram if present; absent columns are simply skipped -- see rerun_msmutect_from_file.py)
LOCUS_COLUMNS = ["CHROMOSOME", "START", "END", "PATTERN",
                 "REFERENCE_SEQUENCE", "REFERENCE_REPEATS", "KNOWN_GERMLINE_VARIANT"]

# locus key used everywhere to compare/join mutations between old and new runs
LOCUS_KEY_COLUMNS = ["CHROMOSOME", "START", "END", "PATTERN"]

# refuse to start a new download/decompression step if free disk space drops below this
MIN_FREE_DISK_GB = 15


def filt_gz_local_path(sample: str) -> str:
    return os.path.join(FILT_DIR, f"{sample}.called.filt.mut.tsv.gz")


def full_gz_local_path(sample: str) -> str:
    return os.path.join(FULL_DIR, f"{sample}.full.mut.tsv.gz")


def rerun_sample_dir(sample: str) -> str:
    return os.path.join(RERUN_DIR, sample)


def rerun_output_prefix(sample: str) -> str:
    return os.path.join(rerun_sample_dir(sample), sample)


def old_yossi_filtered_path(sample: str) -> str:
    return os.path.join(OLD_YOSSI_DIR, f"{sample}.old_yossi_filtered.tsv")


def new_yossi_filtered_path(sample: str) -> str:
    return os.path.join(NEW_YOSSI_DIR, f"{sample}.new_yossi_filtered.tsv")


def ensure_free_disk_gb(min_gb: float = MIN_FREE_DISK_GB, path: str = THIS_DIR):
    import shutil
    free_gb = shutil.disk_usage(path).free / (1024 ** 3)
    if free_gb < min_gb:
        raise RuntimeError(
            f"Only {free_gb:.1f} GB free at {path} (minimum required: {min_gb} GB). "
            f"Aborting to avoid filling up the disk -- clean up data/ or results/ and retry."
        )
