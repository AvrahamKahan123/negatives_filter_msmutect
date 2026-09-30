import csv
import os
import sys

NORMAL_COLS = [f"NORMAL_SUPPORTING_READS_{i}" for i in range(1, 4)]
TUMOR_COLS = [f"TUMOR_SUPPORTING_READS_{i}" for i in range(1, 5)]
MAX_PATTERN_LEN = 15


def to_int(value) -> int:
    """Parse an int, treating NA / empty / junk as 0."""
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 0


def parse_coord(value):
    """Parse a coordinate strictly: None if missing or non-numeric."""
    try:
        return int(str(value).strip())
    except (TypeError, ValueError):
        return None


def total(row: dict, cols: list) -> int:
    return sum(to_int(row.get(c)) for c in cols)


def main(tsv_fp: str):
    callable_loci = [0] * MAX_PATTERN_LEN
    callable_bases = [0] * MAX_PATTERN_LEN
    mutated_loci = [0] * MAX_PATTERN_LEN
    mutated_bases = [0] * MAX_PATTERN_LEN
    bad_coords = 0

    with open(tsv_fp, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if total(row, NORMAL_COLS) < 10 or total(row, TUMOR_COLS) < 10:
                continue

            pattern_len = len(row["PATTERN"].strip())
            if not 1 <= pattern_len <= MAX_PATTERN_LEN:
                continue

            start = parse_coord(row.get("START"))
            end = parse_coord(row.get("END"))+1
            if start is None or end is None or end < start:
                bad_coords += 1
                continue
            length = end - start

            i = pattern_len - 1
            callable_loci[i] += 1
            callable_bases[i] += length
            if row.get("CALL", "").strip() == "M":
                mutated_loci[i] += 1
                mutated_bases[i] += length

    directory, basename = os.path.split(tsv_fp)
    stem = basename
    for suffix in (".gz", ".tsv", ".mut", ".full"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
    out_fp = os.path.join(directory, stem + ".count.txt")

    with open(out_fp, "w") as f:
        f.write("PATTERN_LEN\tCALLABLE_LOCI\tCALLABLE_BASES\t"
                "MUTATED_LOCI\tMUTATED_BASES\n")
        for length_i in range(MAX_PATTERN_LEN):
            f.write(f"{length_i + 1}\t{callable_loci[length_i]}\t"
                    f"{callable_bases[length_i]}\t{mutated_loci[length_i]}\t"
                    f"{mutated_bases[length_i]}\n")

    if bad_coords:
        print(f"warning: skipped {bad_coords} row(s) with missing or invalid "
              f"START/END", file=sys.stderr)
    print(f"Wrote {out_fp}")


if __name__ == '__main__':
    if len(sys.argv) == 2:
        main(sys.argv[1])
    else:
        main("/home/avraham/MaruvkaLab/msmutect_postprocessing/results/GIAB_results_per_pattern_length/hg001_normal1_tumor1.filtered.full.mut.tsv.tsv")