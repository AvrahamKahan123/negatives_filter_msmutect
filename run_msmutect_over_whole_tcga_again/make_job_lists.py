"""Build the input files that the condor submit files `queue ... from`.

Each list is one sample name per line. By default a sample is left OUT if its output for
that step already exists, so re-running this after a partial batch gives you exactly the
work that remains -- submit, wait, re-make, re-submit until the list is empty.

    python -m run_msmutect_over_whole_tcga_again.make_job_lists --step 1
    condor_submit step1_filter.sub

    python -m run_msmutect_over_whole_tcga_again.make_job_lists --step all   # all three

WHERE THE SAMPLE NAMES COME FROM
--------------------------------
Each step enumerates the directory it READS, not the original input_dir:

    step 1  <- input_dir                 (the only step for which input_dir is authoritative)
    step 2  <- <work_root>/filtered
    step 3  <- <work_root>/filtered_wgermline

This matters because samples can enter the pipeline sideways. A sample whose original is
unreadable on the source storage gets recovered to unread/ and filtered by filter_unread.py
-- it is perfectly real from step 2 onwards, but it is NOT in input_dir's listing, so
enumerating input_dir for every step would silently drop it from the rest of the run.

Override that when you need to:

    --from-file PATH   one sample name (or filename) per line
    --from-dir  PATH   every sample in this directory
    --samples   A B C  these samples, named outright
    --all-steps-from-input-dir   the old behaviour: input_dir drives all three steps, and
                                 samples whose input is not ready yet are reported as such
"""

import argparse
import os

from run_msmutect_over_whole_tcga_again import config

STEPS = {
    1: ("step1_filter.txt", "filter"),
    2: ("step2_germline.txt", "germline"),
    3: ("step3_msmutect.txt", "msmutect"),
}

# never mistake a half-written or stray file for a finished sample
IGNORED_SUFFIXES = (".partial", "_unreadable.txt", ".log", ".out", ".err")


def source_dir(step: int) -> str:
    """The directory this step READS -- and therefore the one that defines its work list."""
    if step == 1:
        return config.INPUT_DIR
    if step == 2:
        return config.FILTERED_DIR
    return config.FILTERED_WGERMLINE_DIR


def samples_in_dir(directory: str, suffix: str = None) -> list:
    """Every sample with a file in `directory`, by name.

    os.listdir rather than glob: glob stats each entry to classify it, and on the troubled
    source storage a stat can raise EIO. Reading names out of the directory needs no stat,
    so a sample whose *contents* are unreadable still gets enumerated here -- and is then
    reported properly by write_list's path_state check rather than vanishing.
    """
    if not os.path.isdir(directory):
        return []
    names = set()
    for entry in os.listdir(directory):
        if entry.startswith(".") or entry.endswith(IGNORED_SUFFIXES):
            continue
        if suffix and not entry.endswith(suffix):
            continue
        names.add(config.sample_name(entry))
    return sorted(names)


def samples_in_file(path: str) -> list:
    """Sample names from a one-per-line file. Blank lines and '#' comments are skipped.

    Entries may be bare names or filenames (TCGA-ZB-A963 or TCGA-ZB-A963.full.mut.tsv.gz),
    because these lists get pasted around between the job lists and the fetch lists, which
    use the two different conventions.
    """
    if not os.path.exists(path):
        raise SystemExit(f"ERROR: --from-file not found: {path}")
    names = []
    seen = set()
    with open(path) as f:
        for line in f:
            # tolerate "name<TAB>reason", the shape find_unreadable_files.py writes
            entry = line.split("\t")[0].strip()
            if not entry or entry.startswith("#"):
                continue
            sample = config.sample_name(entry)
            if sample not in seen:
                seen.add(sample)
                names.append(sample)
    if not names:
        raise SystemExit(f"ERROR: no sample names in {path}")
    return names


def all_samples() -> list:
    """Every sample in the source directory, from its *.full.mut.tsv.gz files."""
    return samples_in_dir(config.INPUT_DIR, suffix=config.INPUT_SUFFIX)


def _done(step: int, sample: str) -> bool:
    if step == 1:
        return os.path.exists(config.filtered_path(sample))
    if step == 2:
        return os.path.exists(config.filtered_wgermline_path(sample))
    return os.path.exists(config.msmutect_output_path(sample))


def _input_path(step: int, sample: str) -> str:
    """The file this step reads."""
    if step == 1:
        return config.source_path(sample)
    if step == 2:
        return config.filtered_path(sample)
    return config.filtered_wgermline_path(sample)


def write_list(step: int, samples: list, include_done: bool = False) -> str:
    filename, _ = STEPS[step]
    os.makedirs(config.JOBS_DIR, exist_ok=True)
    # condor will not create these itself, and a missing log dir kills the submit
    os.makedirs(os.path.join(config.LOGS_DIR, f"step{step}"), exist_ok=True)

    ready, missing, unreadable = [], [], []
    for sample in samples:
        state, detail = config.path_state(_input_path(step, sample))
        if state == config.OK:
            ready.append(sample)
        elif state == config.MISSING:
            missing.append(sample)
        else:
            unreadable.append((sample, detail))

    todo = ready if include_done else [s for s in ready if not _done(step, s)]
    done = len(ready) - len(todo)

    path = os.path.join(config.JOBS_DIR, filename)
    with open(path, "w") as f:
        for sample in todo:
            f.write(sample + "\n")

    print(f"step {step}: {len(todo):,} job(s) -> {path}")
    if done:
        print(f"         {done:,} already done (--include-done to redo)")
    if missing:
        where = "input_dir" if step == 1 else f"step {step - 1} has not produced their input"
        print(f"         {len(missing):,} not ready ({where})")

    if unreadable:
        # NOT the same as missing: the file is there but the filesystem cannot read it.
        # Left out of the job list, but recorded -- silently dropping these would give an
        # incomplete final dataset with nothing to show for it.
        bad_path = os.path.join(config.JOBS_DIR, filename.replace(".txt", "_unreadable.txt"))
        with open(bad_path, "w") as f:
            for sample, detail in unreadable:
                f.write(f"{sample}\t{detail}\n")
        print(f"\n  *** {len(unreadable):,} input file(s) EXIST BUT CANNOT BE READ ***")
        for sample, detail in unreadable[:5]:
            print(f"        {sample}  ({detail})")
        if len(unreadable) > 5:
            print(f"        ... and {len(unreadable) - 5:,} more")
        print(f"      full list -> {bad_path}")
        print( "      This is a storage fault, not a pipeline state. These samples are NOT")
        print( "      in the job list and will be missing from the final results until the")
        print(f"      filesystem serves them again; re-run this to pick them up.\n")
    return path


def samples_for_step(step: int, args) -> tuple:
    """-> (samples, description of where they came from)."""
    if args.samples:
        return list(args.samples), "--samples"
    if args.from_file:
        return samples_in_file(args.from_file), args.from_file
    if args.from_dir:
        return samples_in_dir(args.from_dir), args.from_dir
    if args.all_steps_from_input_dir:
        return all_samples(), config.INPUT_DIR
    directory = source_dir(step)
    suffix = config.INPUT_SUFFIX if step == 1 else None
    return samples_in_dir(directory, suffix=suffix), directory


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--step", required=True, choices=["1", "2", "3", "all"])
    parser.add_argument("--include-done", action="store_true",
                        help="list every sample, including ones already finished")
    source = parser.add_argument_group("where the sample names come from "
                                       "(default: the directory the step reads)")
    source.add_argument("--samples", nargs="+", default=None,
                        help="restrict to these samples, named outright")
    source.add_argument("--from-file", metavar="PATH", default=None,
                        help="one sample name or filename per line")
    source.add_argument("--from-dir", metavar="PATH", default=None,
                        help="every sample with a file in this directory")
    source.add_argument("--all-steps-from-input-dir", action="store_true",
                        help="old behaviour: drive all three steps from input_dir, reporting "
                             "samples whose input is not ready yet")
    args = parser.parse_args()

    chosen = [bool(args.samples), bool(args.from_file), bool(args.from_dir),
              args.all_steps_from_input_dir]
    if sum(chosen) > 1:
        parser.error("give at most one of --samples / --from-file / --from-dir / "
                     "--all-steps-from-input-dir")

    steps = [1, 2, 3] if args.step == "all" else [int(args.step)]
    for step in steps:
        samples, where = samples_for_step(step, args)
        print(f"{len(samples):,} sample(s) from {where}")
        if not samples:
            # an empty preceding directory is the normal state before that step has run, so
            # say which directory was empty rather than failing with a generic message
            print(f"step {step}: nothing to do -- {where} has no samples yet\n")
            continue
        write_list(step, samples, include_done=args.include_done)
        print()


if __name__ == "__main__":
    main()
