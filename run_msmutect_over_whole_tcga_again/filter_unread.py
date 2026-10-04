"""Filter the recovered files in unread/ into filtered/, joining them back into the run.

The files in <work_root>/unread/ were fetched from the backup host because the main storage
could not read the originals. They are the same format, so they go through exactly the same
step-1 filter (CALL=="M") and land in the same place — after this, those samples are
indistinguishable from the ones that filtered normally, and steps 2 and 3 pick them up with
no special handling.

    <work_root>/unread/TCGA-ZM-AA0N.full.mut.tsv.gz
        -> <work_root>/filtered/TCGA-ZM-AA0N.full.mut.tsv.gz

Run it over everything (serial, simplest):

    python -m run_msmutect_over_whole_tcga_again.filter_unread

Or one sample, so condor can fan it out (see filter_unread.sub):

    python -m run_msmutect_over_whole_tcga_again.filter_unread --sample TCGA-ZM-AA0N

Already-filtered samples are skipped, so this is safe to re-run as more files are recovered.
A file that is itself damaged is reported and skipped rather than stopping the batch.
"""

import argparse
import glob
import os
import sys

from run_msmutect_over_whole_tcga_again import config, filter_to_mutations

SUFFIXES = ("", ".full.mut.tsv.gz", ".full.mut.tsv")


def recovered_files() -> list:
    """-> [(sample, path)] for everything currently in unread/."""
    found = []
    for suffix in SUFFIXES:
        found.extend(glob.glob(os.path.join(config.UNREAD_DIR, "*" + suffix)))
    by_sample = {}
    for path in sorted(found):
        # a sample present as both .gz and plain would otherwise be filtered twice
        by_sample.setdefault(config.sample_name(path), path)
    return sorted(by_sample.items())


def source_for(sample: str) -> tuple:
    """-> (sample_name, path). Accepts either a bare sample name or a filename.

    jobs/filter_unread.txt holds SAMPLE NAMES while fetch_list.txt holds FILENAMES, so
    normalise rather than fail confusingly when the two get crossed.
    """
    sample = config.sample_name(sample)
    for suffix in SUFFIXES:
        candidate = os.path.join(config.UNREAD_DIR, sample + suffix)
        if os.path.exists(candidate):
            return sample, candidate
    raise SystemExit(f"ERROR: no recovered file for {sample} in {config.UNREAD_DIR}")


def filter_one(sample: str, source: str, force: bool = False, quiet: bool = False) -> str:
    state, detail = config.path_state(source)
    if state == config.UNREADABLE:
        raise RuntimeError(f"recovered copy is itself unreadable ({detail})")
    if state == config.MISSING:
        raise RuntimeError(f"{source} disappeared")
    return filter_to_mutations.filter_sample(sample, force=force, quiet=quiet, source=source)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sample", default=None,
                        help="filter just this one (default: everything in unread/)")
    parser.add_argument("--force", action="store_true",
                        help="re-filter even if the output already exists")
    parser.add_argument("--list", metavar="PATH", default=None,
                        help="write a condor job list of the samples needing work, and exit")
    args = parser.parse_args()

    if not os.path.isdir(config.UNREAD_DIR):
        raise SystemExit(f"ERROR: {config.UNREAD_DIR} does not exist -- nothing recovered yet "
                         f"(see make_fetch_list.sh / fetch_unread.sub)")

    if args.sample:
        samples = [source_for(args.sample)]
    else:
        samples = recovered_files()
        if not samples:
            raise SystemExit(f"no recovered files in {config.UNREAD_DIR}")

    if args.list:
        todo = [s for s, _ in samples
                if args.force or not os.path.exists(config.filtered_path(s))]
        os.makedirs(os.path.dirname(os.path.abspath(args.list)), exist_ok=True)
        os.makedirs(os.path.join(config.LOGS_DIR, "filter_unread"), exist_ok=True)
        with open(args.list, "w") as f:
            for sample in todo:
                f.write(sample + "\n")
        print(f"{len(todo):,} of {len(samples):,} recovered sample(s) need filtering "
              f"-> {args.list}")
        return

    print(f"{len(samples):,} recovered file(s) in {config.UNREAD_DIR}")
    done = skipped = failed = 0
    for i, (sample, source) in enumerate(samples, 1):
        prefix = f"[{i}/{len(samples)}]"
        if os.path.exists(config.filtered_path(sample)) and not args.force:
            print(f"{prefix} {sample}: already filtered, skipping")
            skipped += 1
            continue
        try:
            filter_one(sample, source, force=args.force, quiet=True)
            print(f"{prefix} {sample}: filtered -> {config.filtered_path(sample)}", flush=True)
            done += 1
        except SystemExit as e:
            # filter_sample exits on a bad input; keep going through the rest of the batch
            print(f"{prefix} {sample}: FAILED -- {e}", flush=True)
            failed += 1
        except Exception as e:  # noqa: BLE001 - one damaged file must not stop the batch
            print(f"{prefix} {sample}: FAILED -- {type(e).__name__}: {e}", flush=True)
            failed += 1

    print(f"\n{done:,} filtered, {skipped:,} already done, {failed:,} failed "
          f"-> {config.FILTERED_DIR}")
    if failed:
        print("The failures are files whose recovered copy is also bad; re-fetch them or "
              "raise them with the admin.")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
