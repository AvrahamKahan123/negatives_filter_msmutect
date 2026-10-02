import argparse
import os

import pandas as pd

import config
import msmutect_format
import old_mutation_calls
from results_postprocessing.analyze_full_tcga_mut_set import yossi_filter

# Step 2 (see Plan.md): run the old version's mutation calls through the current
# yossi_filter(), to see which of them the newest post-processing filter would still keep.
#
# The input is the old run's own CALL=="M" rows extracted from its *.full.mut.tsv.gz, NOT
# *.called.filt.mut.tsv.gz: the latter has already been filtered and is missing some of the
# old version's calls (2,442 of them on TCGA-A6-5661), which would bias the comparison in
# step 4. See old_mutation_calls.py.


def run_one(sample: str, noisy_db=None) -> pd.DataFrame:
    src = old_mutation_calls.extract(sample)

    print(f"[{sample}] loading {src}")
    mutations_df = pd.read_csv(src, delimiter="\t", low_memory=False)
    mutations_df = msmutect_format.normalize_for_yossi_filter(mutations_df)
    filtered = yossi_filter(mutations_df, noisy_db or config.noisy_locus_db(), name=sample)

    os.makedirs(config.OLD_YOSSI_DIR, exist_ok=True)
    dest = config.old_yossi_filtered_path(sample)
    filtered.to_csv(dest, sep="\t", index=False)
    print(f"[{sample}] {len(filtered)}/{len(mutations_df)} old calls pass yossi_filter -> {dest}")
    return filtered


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--samples", nargs="+", default=config.SAMPLES)
    args = parser.parse_args()

    for sample in args.samples:
        run_one(sample)


if __name__ == "__main__":
    main()
