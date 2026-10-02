import argparse
import csv
import gzip
import os

import config

# Pulls the rows an MSMuTect run actually called as mutations (CALL == "M") out of a
# *.full.mut.tsv, discarding the ~27M non-mutation loci.
#
# `extract()` is the pipeline's use of this: it caches the OLD run's own mutation calls so
# steps 2 and 4 have a baseline for "what the old version called". That baseline is taken
# from the full file rather than *.called.filt.mut.tsv.gz because the latter has already
# been filtered -- on TCGA-A6-5661 the old full run holds 833,273 M calls but the filt file
# carries only 830,831, and measuring a subset claim against a pre-filtered set would count
# those 2,442 missing calls as new-version violations.
#
# make_filtered_full_old.py uses write_mutation_calls() directly to produce standalone
# filtered copies under data/filtered_full_old/.

MUTATION_CALL = "M"
CALL_COLUMN = "CALL"
PROGRESS_EVERY = 5_000_000


def open_maybe_gzip(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", newline="")
    return open(path, newline="")


def open_write_maybe_gzip(path: str, gzipped: bool):
    # `gzipped` comes from the real destination name, not from `path`, because the file is
    # written to a .partial temporary whose own extension says nothing about the encoding
    if gzipped:
        return gzip.open(path, "wt", newline="\n")
    return open(path, "w", newline="\n")


def write_mutation_calls(source: str, destination: str, progress: bool = True) -> tuple:
    """Stream `source`, writing the header plus every CALL=="M" row to `destination`.

    Returns (kept, total). Writes to a .partial file and renames only on success, so an
    interrupted run never leaves a truncated file that a later cache check would reuse.
    """
    os.makedirs(os.path.dirname(os.path.abspath(destination)), exist_ok=True)
    partial = destination + ".partial"
    gzipped = destination.endswith(".gz")
    kept = total = 0
    with open_maybe_gzip(source) as src, open_write_maybe_gzip(partial, gzipped) as dst:
        reader = csv.reader(src, delimiter="\t")
        writer = csv.writer(dst, delimiter="\t", lineterminator="\n")
        try:
            header = next(reader)
        except StopIteration:
            os.remove(partial)
            raise RuntimeError(f"{source} is empty")
        if CALL_COLUMN not in header:
            os.remove(partial)
            raise RuntimeError(f"{source} has no {CALL_COLUMN} column; is it a *.full.mut.tsv?")
        call_index = header.index(CALL_COLUMN)
        writer.writerow(header)
        for row in reader:
            total += 1
            if row[call_index] == MUTATION_CALL:
                writer.writerow(row)
                kept += 1
            if progress and total % PROGRESS_EVERY == 0:
                print(f"    {total:,} loci scanned, {kept:,} mutations so far...", flush=True)
    os.replace(partial, destination)
    return kept, total


def extract(sample: str, force: bool = False) -> str:
    """Cache the old run's mutation calls for `sample`; returns the cached file's path."""
    destination = config.old_full_mutations_path(sample)
    if os.path.exists(destination) and not force:
        print(f"[{sample}] using cached old-run mutation calls: {destination}")
        return destination

    source = config.full_gz_local_path(sample)
    if not os.path.exists(source):
        raise FileNotFoundError(f"{source} not found -- run download_all.py first")

    print(f"[{sample}] extracting {CALL_COLUMN}=={MUTATION_CALL} rows from {os.path.basename(source)}")
    kept, total = write_mutation_calls(source, destination)
    print(f"[{sample}] {kept:,} mutation calls out of {total:,} loci -> {destination}")
    return destination


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--samples", nargs="+", default=config.SAMPLES)
    parser.add_argument("--force", action="store_true", help="re-extract even if cached")
    args = parser.parse_args()
    for sample in args.samples:
        extract(sample, force=args.force)


if __name__ == "__main__":
    main()
