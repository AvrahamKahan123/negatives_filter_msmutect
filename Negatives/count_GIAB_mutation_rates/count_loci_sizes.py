"""Count rows by PATTERN length across every TSV in a directory.

Reads each *.tsv (and *.tsv.gz) in the given directory, tallies rows by the
length of the PATTERN field, and adds the per-file tallies into one total.

    python3 count_pattern_lengths_dir.py /path/to/dir
    python3 count_pattern_lengths_dir.py /path/to/dir -o totals.tsv --per-file
"""
import argparse
import csv
import gzip
import os
import sys
from collections import Counter

PATTERNS = ("*.tsv", "*.tsv.gz")


def open_maybe_gz(fp: str):
    if fp.endswith(".gz"):
        return gzip.open(fp, "rt", newline="")
    return open(fp, newline="")


def count_file(fp: str, col: str) -> Counter:
    """Tally rows in one file by the length of `col`."""
    counts = Counter()
    with open_maybe_gz(fp) as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("file is empty")
        fields = [name.strip() for name in reader.fieldnames]
        if col not in fields:
            raise ValueError(f"no column named {col!r} (found: {', '.join(fields[:6])})")
        reader.fieldnames = fields
        for row in reader:
            value = row.get(col) or ""
            counts[len(value.strip())] += 1
    return counts


def find_inputs(directory: str) -> list:
    import glob
    found = []
    for pattern in PATTERNS:
        found.extend(glob.glob(os.path.join(directory, pattern)))
    return sorted(set(found))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", help="directory containing the TSV files")
    parser.add_argument("-o", "--output", help="write totals here (default: stdout)")
    parser.add_argument("-c", "--col", default="PATTERN", help="column name (default: PATTERN)")
    parser.add_argument("--per-file", action="store_true",
                        help="also print a per-file breakdown to stderr")
    args = parser.parse_args()

    inputs = find_inputs(args.directory)
    if not inputs:
        sys.exit(f"no .tsv or .tsv.gz files found in {args.directory}")

    totals = Counter()
    n_ok = 0
    for fp in inputs:
        try:
            counts = count_file(fp, args.col)
        except (OSError, ValueError) as err:
            print(f"skipping {os.path.basename(fp)}: {err}", file=sys.stderr)
            continue
        totals += counts
        n_ok += 1
        if args.per_file:
            rows = sum(counts.values())
            print(f"{os.path.basename(fp)}\t{rows} rows\t"
                  f"lengths {min(counts)}-{max(counts)}" if counts else
                  f"{os.path.basename(fp)}\t0 rows", file=sys.stderr)

    if not totals:
        sys.exit("no rows counted")

    out = open(args.output, "w") if args.output else sys.stdout
    try:
        out.write("PATTERN_LEN\tN_ROWS\n")
        for length in range(0, max(totals) + 1):
            if length in totals:
                out.write(f"{length}\t{totals[length]}\n")
    finally:
        if args.output:
            out.close()

    print(f"Counted {sum(totals.values())} rows from {n_ok}/{len(inputs)} file(s)",
          file=sys.stderr)


if __name__ == "__main__":
    main()