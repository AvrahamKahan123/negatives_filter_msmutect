import argparse

import config
import download_samples
import run_yossi_filter_on_old
from Negatives.run_all_samples_through_newest_filter import rerun_msmutect_from_file
import compare_old_vs_new
import ssh_connection

# Runs the pipeline described in Plan.md over one or more samples, in this process:
#   1. make sure the sample's filt + full files are present (download if not)
#   2. run the old filt file through yossi_filter()
#   3. rerun the newest MSMuTect on the full file via --from_file
#   4. compare the old vs. new (raw + yossi_filter'd) mutation sets
#
# This is also the per-sample unit of work that run_parallel.py spawns -- it invokes
#     python run_all.py --samples <ONE_SAMPLE> --skip-summary
# once per sample, concurrently. Each sample writes only its own files (step 4 saves a
# per-sample row via compare_old_vs_new.save_row), so concurrent workers never collide;
# run_parallel.py builds summary.tsv from those rows once they have all finished.
#
# Run `conda activate genomics` first (pandas/numpy/pysam/scipy all live there).
#
# Prefer running download_all.py before this: with everything already local, the pipeline
# needs no network at all and cannot stall on an expired SSH connection.


def run_sample(sample: str, force: bool = False) -> dict:
    print(f"\n===== {sample} =====")
    download_samples.download_filt(sample, force=force)
    run_yossi_filter_on_old.run_one(sample)
    rerun_msmutect_from_file.run_sample(sample, force=force)
    row = compare_old_vs_new.compare_one(sample)
    compare_old_vs_new.save_row(row)
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--samples", nargs="+", default=config.SAMPLES)
    parser.add_argument("--force", action="store_true", help="redo every step even if outputs already exist")
    parser.add_argument("--skip-summary", action="store_true",
                         help="do not write summary.tsv (set by run_parallel.py, which writes it "
                              "once after every worker has finished)")
    parser.add_argument("--close-ssh", action="store_true",
                         help="close the shared SSH connection when finished")
    args = parser.parse_args()

    for sample in args.samples:
        run_sample(sample, force=args.force)

    if not args.skip_summary:
        print("\n===== summary =====")
        compare_old_vs_new.summarize(args.samples)

    if args.close_ssh:
        ssh_connection.close_master()


if __name__ == "__main__":
    main()
