import csv
import os
import sys
from collections import defaultdict

NORMAL_COLS = [f"NORMAL_SUPPORTING_READS_{i}" for i in range(1, 4)]
TUMOR_COLS = [f"TUMOR_SUPPORTING_READS_{i}" for i in range(1, 5)]
MAX_PATTERN_LEN = 15


def to_int(value) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 0


def total(row: dict, cols: list) -> int:
    return sum(to_int(row.get(c)) for c in cols)


def main(tsv_fp: str):
    # count separately for each (pattern length, number of repeats) combination
    callable_counts = defaultdict(int)
    mutated_counts = defaultdict(int)

    with open(tsv_fp, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if total(row, NORMAL_COLS) < 10 or total(row, TUMOR_COLS) < 10:
                continue
            pattern_len = len(row["PATTERN"].strip())
            if not 1 <= pattern_len <= MAX_PATTERN_LEN:
                continue
            num_repeats = to_int(row["REFERENCE_REPEATS"])
            key = (pattern_len, num_repeats)
            callable_counts[key] += 1
            if row.get("CALL", "").strip() == "M":
                mutated_counts[key] += 1

    directory, basename = os.path.split(tsv_fp)
    out_fp = os.path.join(directory, basename.split(".")[0] + ".count.txt")
    with open(out_fp, "w") as f:
        f.write("PATTERN_LEN\tNUM_REPEATS\tCALLABLE\tMUTATED\n")
        for pattern_len, num_repeats in sorted(callable_counts):
            key = (pattern_len, num_repeats)
            f.write(f"{pattern_len}\t{num_repeats}\t{callable_counts[key]}\t{mutated_counts[key]}\n")

    print(f"Wrote {out_fp}")


if __name__ == '__main__':
    if len(sys.argv) == 2:
        main(sys.argv[1])
    else:
        main("/home/avraham/MaruvkaLab/msmutect_postprocessing/data/gib_files_filtered/hg002_normal2_tumor1.filtered.full.mut.tsv")
