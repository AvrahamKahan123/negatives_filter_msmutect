"""Merge the per-sample .count.txt files produced by count_callable.py.

Each input file has one row per pattern length (1-15) with CALLABLE and
MUTATED columns. This script combines any number of them into a single
table, either one row per (sample, pattern_length) or summed across all
samples.
"""
import argparse
import os
import sys

MAX_PATTERN_LEN = 15
HEADER = ("PATTERN_LEN", "CALLABLE", "MUTATED")


def sample_name(fp: str) -> str:
    """HG001_msmutect_normal0_tumor0.count.txt -> HG001_msmutect_normal0_tumor0"""
    return os.path.basename(fp).split(".")[0]


def read_counts(fp: str):
    """Return {pattern_len: (callable, mutated)}. Tolerates a missing header."""
    counts = {}
    with open(fp) as f:
        for line_no, line in enumerate(f, start=1):
            fields = line.strip().split("\t")
            if len(fields) < 3:
                if line.strip():
                    print(f"{fp}:{line_no}: skipping malformed line", file=sys.stderr)
                continue
            if fields[0].upper() == HEADER[0]:
                continue
            try:
                length, callable_n, mutated_n = (int(x) for x in fields[:3])
            except ValueError:
                print(f"{fp}:{line_no}: skipping non-numeric line", file=sys.stderr)
                continue
            counts[length] = (callable_n, mutated_n)
    if not counts:
        raise ValueError(f"no usable rows found in {fp}")
    return counts


def rate(callable_n: int, mutated_n: int) -> str:
    return f"{mutated_n / callable_n:.6g}" if callable_n else "NA"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("count_files", nargs="+", help="one or more .count.txt files")
    parser.add_argument("-o", "--output", default="merged_counts.tsv",
                        help="output path (default: merged_counts.tsv)")
    parser.add_argument("--aggregate", action="store_true",
                        help="sum across samples instead of one row per sample")
    args = parser.parse_args()

    per_sample = {}
    for fp in args.count_files:
        try:
            per_sample[sample_name(fp)] = read_counts(fp)
        except (OSError, ValueError) as err:
            print(f"skipping {fp}: {err}", file=sys.stderr)

    if not per_sample:
        sys.exit("no readable input files")

    lengths = sorted({n for counts in per_sample.values() for n in counts}
                     | set(range(1, MAX_PATTERN_LEN + 1)))

    with open(args.output, "w") as out:
        if args.aggregate:
            out.write("PATTERN_LEN\tCALLABLE\tMUTATED\tMUTATION_RATE\tN_SAMPLES\n")
            for length in lengths:
                rows = [c[length] for c in per_sample.values() if length in c]
                callable_n = sum(r[0] for r in rows)
                mutated_n = sum(r[1] for r in rows)
                out.write(f"{length}\t{callable_n}\t{mutated_n}\t"
                          f"{rate(callable_n, mutated_n)}\t{len(rows)}\n")
        else:
            out.write("SAMPLE\tPATTERN_LEN\tCALLABLE\tMUTATED\tMUTATION_RATE\n")
            for sample in sorted(per_sample):
                counts = per_sample[sample]
                for length in lengths:
                    callable_n, mutated_n = counts.get(length, (0, 0))
                    out.write(f"{sample}\t{length}\t{callable_n}\t{mutated_n}\t"
                              f"{rate(callable_n, mutated_n)}\n")

    print(f"Merged {len(per_sample)} sample(s) -> {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
