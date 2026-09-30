import argparse
import os

import pandas as pd

import config
from results_postprocessing.NoisyLocusDB import NoisyLocusDB
from results_postprocessing.analyze_full_tcga_mut_set import yossi_filter

# Step 2 (see Plan.md): run every downloaded *.called.filt.mut.tsv.gz (the old version's
# mutation calls) through the current yossi_filter(), to see which of those old calls the
# newest post-processing filter would still keep.


def run_one(sample: str) -> pd.DataFrame:
    src = config.filt_gz_local_path(sample)
    if not os.path.exists(src):
        raise FileNotFoundError(f"{src} not found -- run download_samples.py first")

    print(f"[{sample}] loading {src}")
    mutations_df = pd.read_csv(src, delimiter="\t")
    noisy_db = NoisyLocusDB()
    filtered = yossi_filter(mutations_df, noisy_db, name=sample)

    os.makedirs(config.OLD_YOSSI_DIR, exist_ok=True)
    dest = config.old_yossi_filtered_path(sample)
    filtered.to_csv(dest, sep="\t", index=False)
    print(f"[{sample}] {len(filtered)}/{len(mutations_df)} old calls pass yossi_filter -> {dest}")
    return filtered


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", nargs="+", default=config.SAMPLES)
    args = parser.parse_args()

    for sample in args.samples:
        run_one(sample)


if __name__ == "__main__":
    main()
