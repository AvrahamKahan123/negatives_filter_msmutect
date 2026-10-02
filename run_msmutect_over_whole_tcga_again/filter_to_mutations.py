"""Step 1: keep only the rows the old run called a mutation.

One condor job per sample:

    <input_dir>/TCGA-W5-AA2K.full.mut.tsv.gz  ->  <work_root>/filtered/TCGA-W5-AA2K.full.mut.tsv.gz

The header is kept, so each output is a valid *.full.mut.tsv in its own right. Everything
else -- NM, AN, INS, RR, FFT, TMA, LOH -- is dropped, which is ~97-99% of the ~27.7M rows.

CALL is located BY NAME rather than by position, because the column index shifts between
format versions (step 2 drops two columns and adds one).

    python -m run_msmutect_over_whole_tcga_again.filter_to_mutations --sample TCGA-W5-AA2K
"""

import argparse
import csv
import gzip
import os
import sys

from run_msmutect_over_whole_tcga_again import config

MUTATION_CALL = "M"
CALL_COLUMN = "CALL"
PROGRESS_EVERY = 5_000_000


def _open_read(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", newline="")
    return open(path, "r", newline="")


def _open_write(path: str, gzipped: bool):
    # `gzipped` comes from the real destination name, not `path`, because output is written
    # to a .partial whose own extension says nothing about the encoding
    if gzipped:
        return gzip.open(path, "wt", newline="")
    return open(path, "w", newline="")


def filter_sample(sample: str, force: bool = False, quiet: bool = False) -> str:
    source = config.source_path(sample)
    destination = config.filtered_path(sample)
    state, detail = config.path_state(source)
    if state == config.UNREADABLE:
        # distinguish a storage fault from an absent file -- os.path.exists() calls both False
        raise SystemExit(f"ERROR: input exists but cannot be read: {source}\n"
                         f"       {detail}\n"
                         f"       Storage fault, not a pipeline problem. Re-submit this sample "
                         f"once the filesystem serves it again.")
    if state == config.MISSING:
        raise SystemExit(f"ERROR: input not found: {source}")
    if os.path.exists(destination) and not force:
        if not quiet:
            print(f"[{sample}] already filtered, skipping (--force to redo): {destination}")
        return destination

    os.makedirs(os.path.dirname(destination), exist_ok=True)
    # write to .partial and rename only on success, so an interrupted or evicted condor job
    # never leaves a truncated file that the skip-existing check would happily reuse
    partial = destination + ".partial"
    kept = total = 0
    try:
        with _open_read(source) as src, _open_write(partial, destination.endswith(".gz")) as dst:
            reader = csv.reader(src, delimiter="\t")
            writer = csv.writer(dst, delimiter="\t", lineterminator="\n")
            try:
                header = next(reader)
            except StopIteration:
                raise SystemExit(f"ERROR: {source} is empty")
            if CALL_COLUMN not in header:
                raise SystemExit(f"ERROR: {source} has no {CALL_COLUMN} column; is it really a "
                                 f"*.full.mut.tsv? (columns: {', '.join(header[:8])}...)")
            call_index = header.index(CALL_COLUMN)
            writer.writerow(header)
            for row in reader:
                total += 1
                if row[call_index] == MUTATION_CALL:
                    writer.writerow(row)
                    kept += 1
                if not quiet and total % PROGRESS_EVERY == 0:
                    print(f"  {total:,} loci scanned, {kept:,} mutations", flush=True)
    except BaseException:
        if os.path.exists(partial):
            os.remove(partial)
        raise

    os.replace(partial, destination)
    if not quiet:
        percent = 100 * kept / max(total, 1)
        print(f"[{sample}] {kept:,} mutations of {total:,} loci ({percent:.2f}%) -> {destination}")
    return destination


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sample", required=True, help="e.g. TCGA-W5-AA2K")
    parser.add_argument("--force", action="store_true", help="re-filter even if the output exists")
    args = parser.parse_args()
    filter_sample(args.sample, force=args.force)


if __name__ == "__main__":
    main()
