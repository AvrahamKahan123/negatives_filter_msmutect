"""Count rows with usable NORMAL_SUPPORTING_READS, split by PATTERN length.

Walks a TSV and keeps two counters: one for rows whose PATTERN is a single
character, one for rows whose PATTERN is two characters. A row is counted
only if its normal supporting reads value is present, numeric, and nonzero.
Every other row is ignored.

    python3 count_supported.py file.tsv
"""
import argparse
import csv
import gzip
import sys

# Some files carry a single NORMAL_SUPPORTING_READS column; others (like the
# msmutect .full.mut.tsv output) split it across numbered columns. Support both.
SINGLE_COL = "NORMAL_SUPPORTING_READS"
NUMBERED_COLS = [f"TUMOR_SUPPORTING_READS_{i}" for i in range(2, 4)]


def open_maybe_gz(fp: str):
    if fp.endswith(".gz"):
        return gzip.open(fp, "rt", newline="")
    return open(fp, newline="")


def to_int(value) -> int:
    """Parse an int, treating NA / empty / junk as 0."""
    try:
        return int(str(value).strip())
    except (TypeError, ValueError):
        return 0


def pick_columns(fieldnames: list) -> list:
    """Return the normal-read columns actually present in this file."""
    if SINGLE_COL in fieldnames:
        print("XXXXXX XXXXXX XXXXXX")
        exit()
        return [SINGLE_COL]
    present = [c for c in NUMBERED_COLS if c in fieldnames]
    if present:
        return present
    sys.exit(f"no {SINGLE_COL} or {SINGLE_COL}_N columns found; "
             f"header is: {', '.join(fieldnames[:8])}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tsv_fp")
    args = parser.parse_args()

    counted = {1: 0, 2: 0}     # rows with usable normal reads
    examined = {1: 0, 2: 0}    # rows of that pattern length, before the read check

    with open_maybe_gz(args.tsv_fp) as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames is None:
            sys.exit(f"{args.tsv_fp} is empty")
        reader.fieldnames = [name.strip() for name in reader.fieldnames]
        if "PATTERN" not in reader.fieldnames:
            sys.exit(f"no PATTERN column; header is: {', '.join(reader.fieldnames[:8])}")
        read_cols = pick_columns(reader.fieldnames)

        for row in reader:
            pattern_len = len((row.get("PATTERN") or "").strip())
            if pattern_len not in counted:
                continue
            examined[pattern_len] += 1
            if sum(to_int(row.get(c)) for c in read_cols) > 0:
                counted[pattern_len] += 1

    print(f"reads column(s): {', '.join(read_cols)}", file=sys.stderr)
    print("PATTERN_LEN\tWITH_READS\tTOTAL_ROWS")
    for length in (1, 2):
        print(f"{length}\t{counted[length]}\t{examined[length]}")


if __name__ == "__main__":
    main()