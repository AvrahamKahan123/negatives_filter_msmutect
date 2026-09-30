import os
import re

import pandas as pd

# all_heterozygous_loci.csv columns:
#   (index), CHROMOSOME, START, END, PATTERN, REFERENCE_SEQUENCE, REFERENCE_REPEATS, ALLELES_REPEATS, HETEROZYGOUS
# ALLELES_REPEATS maps a sample's normal genotype (allele_1, allele_2) -> number of samples.
#   (a, a) homozygous ;  (a, b) with a != b heterozygous ;  (0, 0) non-call.
# The existing HETEROZYGOUS flag also counts loci whose only non-reference signal is
# non-calls, so here we recompute heterozygosity directly from the genotypes.

DATA_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT = os.path.join(DATA_DIR, "all_heterozygous_loci.csv")
OUT_ACTUAL = os.path.join(DATA_DIR, "actual_heterozygous_loci.csv")
OUT_MIN_SAMPLES = os.path.join(DATA_DIR, "actual_heterozygous_loci_min200_samples.csv")

CHUNKSIZE = 500_000
MIN_SAMPLES = 200

GENOTYPE_RE = re.compile(r"\((\d+),\s*(\d+)\):\s*(\d+)")


def heterozygous_sample_count(alleles_repeats) -> int:
    # number of samples that genuinely support multiple versions (two distinct, non-zero alleles),
    # i.e. ignoring homozygous (a, a) and non-call (0, 0) genotypes.
    if not isinstance(alleles_repeats, str):
        return 0
    total = 0
    for a, b, count in GENOTYPE_RE.findall(alleles_repeats):
        a, b = int(a), int(b)
        if a != b and a != 0 and b != 0:
            total += int(count)
    return total


def main():
    first_chunk = True
    n_total = n_actual = n_min_samples = 0

    for chunk in pd.read_csv(INPUT, index_col=0, chunksize=CHUNKSIZE):
        chunk["HET_SAMPLE_COUNT"] = chunk["ALLELES_REPEATS"].map(heterozygous_sample_count)

        actually_het = chunk[chunk["HET_SAMPLE_COUNT"] > 0]           # at least one heterozygous sample
        min_samples = chunk[chunk["HET_SAMPLE_COUNT"] >= MIN_SAMPLES]  # >= 200 heterozygous samples

        actually_het.to_csv(OUT_ACTUAL, mode="w" if first_chunk else "a", header=first_chunk, index=False)
        min_samples.to_csv(OUT_MIN_SAMPLES, mode="w" if first_chunk else "a", header=first_chunk, index=False)

        n_total += len(chunk)
        n_actual += len(actually_het)
        n_min_samples += len(min_samples)
        first_chunk = False
        print(f"processed {n_total:,} loci -> actually het {n_actual:,} -> >={MIN_SAMPLES} het samples {n_min_samples:,}")

    print(f"\nDONE\n  input loci:              {n_total:,}"
          f"\n  actually heterozygous:   {n_actual:,}  -> {OUT_ACTUAL}"
          f"\n  >= {MIN_SAMPLES} het samples:     {n_min_samples:,}  -> {OUT_MIN_SAMPLES}")


if __name__ == "__main__":
    main()
