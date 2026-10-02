"""Build the inputs files that the condor submit files `queue ... from`.

Each list is one sample name per line. By default a sample is left OUT if its output for
that step already exists, so re-running this after a partial batch gives you exactly the
work that remains -- submit, wait, re-make, re-submit until the list is empty.

    python -m run_msmutect_over_whole_tcga_again.make_job_lists --step 1
    condor_submit step1_filter.sub

    python -m run_msmutect_over_whole_tcga_again.make_job_lists --step all   # all three
"""

import argparse
import glob
import os

from run_msmutect_over_whole_tcga_again import config

STEPS = {
    1: ("step1_filter.txt", "filter"),
    2: ("step2_germline.txt", "germline"),
    3: ("step3_msmutect.txt", "msmutect"),
}


def all_samples() -> list:
    """Every sample in the source directory, from its *.full.mut.tsv.gz files."""
    pattern = os.path.join(config.INPUT_DIR, "*" + config.INPUT_SUFFIX)
    return sorted({config.sample_name(p) for p in glob.glob(pattern)})


def _done(step: int, sample: str) -> bool:
    if step == 1:
        return os.path.exists(config.filtered_path(sample))
    if step == 2:
        return os.path.exists(config.filtered_wgermline_path(sample))
    return os.path.exists(config.msmutect_output_path(sample))


def _ready(step: int, sample: str) -> bool:
    """Is this step's INPUT present? Keeps a list from containing jobs that must fail."""
    if step == 1:
        return os.path.exists(config.source_path(sample))
    if step == 2:
        return os.path.exists(config.filtered_path(sample))
    return os.path.exists(config.filtered_wgermline_path(sample))


def write_list(step: int, samples: list, include_done: bool = False) -> str:
    filename, _ = STEPS[step]
    os.makedirs(config.JOBS_DIR, exist_ok=True)
    # condor will not create these itself, and a missing log dir kills the submit
    os.makedirs(os.path.join(config.LOGS_DIR, f"step{step}"), exist_ok=True)

    ready = [s for s in samples if _ready(step, s)]
    not_ready = len(samples) - len(ready)
    todo = ready if include_done else [s for s in ready if not _done(step, s)]
    done = len(ready) - len(todo)

    path = os.path.join(config.JOBS_DIR, filename)
    with open(path, "w") as f:
        for sample in todo:
            f.write(sample + "\n")

    print(f"step {step}: {len(todo):,} job(s) -> {path}")
    if done:
        print(f"         {done:,} already done (--include-done to redo)")
    if not_ready:
        print(f"         {not_ready:,} not ready (step {step - 1} has not produced their input)")
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--step", required=True, choices=["1", "2", "3", "all"])
    parser.add_argument("--include-done", action="store_true",
                        help="list every sample, including ones already finished")
    parser.add_argument("--samples", nargs="+", default=None,
                        help="restrict to these samples (default: everything in input_dir)")
    args = parser.parse_args()

    samples = args.samples or all_samples()
    if not samples:
        raise SystemExit(f"no samples found: no *{config.INPUT_SUFFIX} in {config.INPUT_DIR}")
    print(f"{len(samples):,} sample(s) in {config.INPUT_DIR}\n")

    steps = [1, 2, 3] if args.step == "all" else [int(args.step)]
    for step in steps:
        write_list(step, samples, include_done=args.include_done)


if __name__ == "__main__":
    main()
