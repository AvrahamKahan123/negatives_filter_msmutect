import os
import sys

# NoisyLocusDB lives in the results_postprocessing package, which mixes package and flat
# imports; put both its package root and its own dir on the path so it imports cleanly.
POSTPROC_DIR = "/"
sys.path.insert(0, POSTPROC_DIR)
sys.path.insert(0, os.path.join(POSTPROC_DIR, "results_postprocessing"))

from converting_from_old_format_to_new.NoisyLocusDB import NoisyLocusDB

loci_file = "/home/avraham/MaruvkaLab/msmutect_development/data/GRCh38.d1.vd1_1to15_repetitive_loci_sorted_fixed"
output_file = loci_file + ".with_germline_variance"

# 0-based columns in the loci file, matching how NoisyLocusDB.create_db encoded its keys.
CHROM_COL, START_COL, STOP_COL, PATTERN_COL = 0, 3, 4, 12

HETEROGENEOUS = "HETEROGENEOUS"  # locus is noisy (varies in the germline population)
HOMOGENEOUS = "HOMOGENEOUS"      # locus is not noisy


def main():
    noisy_db = NoisyLocusDB()  # O(1) set-backed membership lookups
    n_total = n_heterogeneous = 0

    with open(loci_file, "r") as infile, open(output_file, "w") as outfile:
        for line in infile:
            line = line.rstrip("\n")
            fields = line.split("\t")
            chromosome = fields[CHROM_COL].replace("chr", "")
            is_noisy = noisy_db.locus_present_in_db(
                chromosome, fields[START_COL], fields[STOP_COL], fields[PATTERN_COL]
            )
            label = HETEROGENEOUS if is_noisy else HOMOGENEOUS
            outfile.write(f"{line}\t{label}\n")
            n_total += 1
            n_heterogeneous += is_noisy

    print(f"wrote {output_file}"
          f"\n  total loci:    {n_total:,}"
          f"\n  HETEROGENEOUS: {n_heterogeneous:,}"
          f"\n  HOMOGENEOUS:   {n_total - n_heterogeneous:,}")


if __name__ == "__main__":
    main()
