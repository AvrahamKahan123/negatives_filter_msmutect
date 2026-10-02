"""Step 3: re-call each sample with the current MSMuTect via --from_file.

One condor job per sample:

    <work_root>/filtered_wgermline/TCGA-W5-AA2K.full.mut.tsv.gz
        -> <work_root>/msmutect_results/TCGA-W5-AA2K.full.mut.tsv  (+ .mutated_only, .call_counts)

A thin wrapper over the running_from_file pyramid: layer 2 splits the file into
pseudo-normal/pseudo-tumor histograms (layer 1) and runs msmutect on them. This module only
supplies the cluster paths.

Because the input is the step-1 filtered file, this re-calls exactly the loci the OLD run
called M -- which is what makes the old-vs-new comparison a subset test.

    python -m run_msmutect_over_whole_tcga_again.run_msmutect_from_file --sample TCGA-W5-AA2K
"""

import argparse
import os

from run_msmutect_over_whole_tcga_again import config  # also puts repo_root on sys.path

from running_from_file import run_single_sample
from running_from_file.split_histograms import MissingGermlineColumnError


def run_sample(sample: str, force: bool = False, keep_histograms: bool = False,
               quiet: bool = False) -> str:
    source = config.filtered_wgermline_path(sample)
    output_prefix = config.msmutect_output_prefix(sample)
    produced = config.msmutect_output_path(sample)
    if not os.path.exists(source):
        raise SystemExit(f"ERROR: {source} not found -- has step 2 run for this sample?")
    if os.path.exists(produced) and not force:
        if not quiet:
            print(f"[{sample}] already re-called, skipping (--force to redo): {produced}")
        return produced

    # fail fast and clearly if the cluster msmutect path in config.json is wrong
    config.require(config.MSMUTECT_PATH, "msmutect_path")
    os.makedirs(os.path.dirname(output_prefix), exist_ok=True)

    try:
        out = run_single_sample.run_from_full_mut_tsv(
            source, output_prefix, msmutect=config.MSMUTECT_PATH,
            keep_histograms=keep_histograms, quiet=quiet)
    except MissingGermlineColumnError as e:
        raise SystemExit(f"ERROR: {e}\n       Step 2 (add_germline) has not been run on this "
                         f"sample, or produced an unannotated file.")

    if not quiet:
        counts = run_single_sample.call_counts(output_prefix)
        if counts:
            print(f"[{sample}] call counts: "
                  + ", ".join(f"{call}={count:,}" for call, count in counts.items()))
        print(f"[{sample}] done -> {out}")
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sample", required=True, help="e.g. TCGA-W5-AA2K")
    parser.add_argument("--force", action="store_true", help="re-run even if the output exists")
    parser.add_argument("--keep-histograms", action="store_true",
                        help="keep the split histograms instead of deleting them on success")
    args = parser.parse_args()
    run_sample(args.sample, force=args.force, keep_histograms=args.keep_histograms)


if __name__ == "__main__":
    main()
