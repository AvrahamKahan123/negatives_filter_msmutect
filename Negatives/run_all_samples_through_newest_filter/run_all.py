import argparse
import os

import pandas as pd

import config
import download_samples
import run_yossi_filter_on_old
import rerun_msmutect_from_file
import compare_old_vs_new
import ssh_connection
from results_postprocessing.NoisyLocusDB import NoisyLocusDB

# Runs the full pipeline described in Plan.md, one sample at a time so at most one
# sample's *.full.mut.tsv.gz is ever being downloaded/processed at once:
#   1. download the sample's filt + full files
#   2. run the old filt file through yossi_filter()
#   3. rerun the newest MSMuTect on the full file via --from_file
#   4. compare the old vs. new (raw + yossi_filter'd) mutation sets
#
# Run `conda activate genomics` first (pandas/numpy/pysam/scipy all live there).
#
# Your SSH password is asked for exactly once, on the first transfer: every download
# shares one SSH connection (see ssh_connection.py), which outlives the individual
# rsyncs and the long MSMuTect runs in between them.


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", nargs="+", default=config.SAMPLES)
    parser.add_argument("--force", action="store_true", help="redo every step even if outputs already exist")
    parser.add_argument("--close-ssh", action="store_true",
                         help="close the shared SSH connection when the pipeline finishes, instead of "
                              "leaving it open for a follow-up run (it expires on its own either way)")
    args = parser.parse_args()

    noisy_db = NoisyLocusDB()
    for sample in args.samples:
        print(f"\n===== {sample} =====")
        download_samples.download_filt(sample, force=args.force)
        run_yossi_filter_on_old.run_one(sample)
        rerun_msmutect_from_file.run_sample(sample, force=args.force)

    print("\n===== comparison =====")
    rows = [compare_old_vs_new.compare_one(sample, noisy_db) for sample in args.samples]

    os.makedirs(config.COMPARISON_DIR, exist_ok=True)
    summary_fp = os.path.join(config.COMPARISON_DIR, "summary.tsv")
    pd.DataFrame(rows).to_csv(summary_fp, sep="\t", index=False)
    print(f"\nsummary -> {summary_fp}")

    overall_ok = all(r["RAW_SUBSET_OK"] and r["YOSSI_SUBSET_OK"] for r in rows)
    print("OVERALL: " + ("PASS -- newest MSMuTect's mutations are a subset of the old version's for every sample"
                          if overall_ok else
                          "FAIL -- see *.raw_violations.tsv / *.yossi_violations.tsv for the offending loci"))

    if args.close_ssh:
        ssh_connection.close_master()


if __name__ == "__main__":
    main()
