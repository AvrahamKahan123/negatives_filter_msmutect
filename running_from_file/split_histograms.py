"""Layer 1 of the pyramid: split a *.full.mut.tsv into pseudo-normal/pseudo-tumor histograms.

A *.full.mut.tsv already holds the per-locus normal and tumor histograms, so an existing
MSMuTect run can be re-called without going back to the BAMs: carve the two histograms back
out and feed them to `msmutect --from_file`.

This layer ONLY splits. The input must already carry KNOWN_GERMLINE_VARIANT -- adding that
column is the job of `converting_from_old_format_to_new`, and this module refuses to run
without it rather than inventing the value. Columns are selected BY NAME, so the split is
robust to columns being added, removed or reordered.

Nothing here runs MSMuTect; `run_single_sample.py` does that.
"""

import argparse
import csv
import gzip
import os
import sys

# Locus columns in output order. KNOWN_GERMLINE_VARIANT sits right after REFERENCE_REPEATS,
# matching what current MSMuTect writes.
LOCUS_COLUMNS = ["CHROMOSOME", "START", "END", "PATTERN",
                 "REFERENCE_SEQUENCE", "REFERENCE_REPEATS", "KNOWN_GERMLINE_VARIANT"]
GERMLINE_COLUMN = "KNOWN_GERMLINE_VARIANT"

# the current --from_file run hard-errors without this in the NORMAL histogram
# (MSMuTect src/Entry/PairFileBatches.py::run_from_file)
CONVERTER_HINT = (
    "Add it first with converting_from_old_format_to_new/convert_old_msmutect_file_to_new.py, "
    "which drops the old 'Noisy Locus' columns and inserts KNOWN_GERMLINE_VARIANT after "
    "REFERENCE_REPEATS."
)

PROGRESS_EVERY = 2_000_000

NORMAL, TUMOR = "NORMAL_", "TUMOR_"


class MissingGermlineColumnError(RuntimeError):
    """Raised when the input predates KNOWN_GERMLINE_VARIANT, so it must be converted first."""


def open_maybe_gzip(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", newline="")
    return open(path, newline="")


def histogram_columns(prefix: str):
    return ([f"{prefix}MOTIF_REPEATS_{i}" for i in range(1, 7)] +
            [f"{prefix}SUPPORTING_READS_{i}" for i in range(1, 7)])


def column_plan(header_index: dict, wanted_columns):
    """Source indices for `wanted_columns`, in order, skipping any not present."""
    return [(name, header_index[name]) for name in wanted_columns if name in header_index]


def histogram_name(full_mut_tsv: str, output_dir: str, which: str) -> str:
    """Conventional histogram path for an input file: <name>.pseudo_{normal,tumor}.hist.tsv"""
    return os.path.join(output_dir, f"{sample_name(full_mut_tsv)}.pseudo_{which}.hist.tsv")


def sample_name(fp: str) -> str:
    name = os.path.basename(fp)
    if name.endswith(".gz"):
        name = name[: -len(".gz")]
    for suffix in (".filtered.full.mut.tsv", ".full.mut.tsv"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return os.path.splitext(name)[0]


def split_full_mut_tsv(full_mut_tsv: str, normal_destination: str, tumor_destination: str,
                       progress: bool = True) -> int:
    """Write both histograms in ONE pass over the input. Returns the number of loci written.

    Raises MissingGermlineColumnError if the input has no KNOWN_GERMLINE_VARIANT column.
    """
    if not os.path.exists(full_mut_tsv):
        raise FileNotFoundError(f"input not found: {full_mut_tsv}")
    for destination in (normal_destination, tumor_destination):
        os.makedirs(os.path.dirname(os.path.abspath(destination)), exist_ok=True)

    with open_maybe_gzip(full_mut_tsv) as source, \
            open(normal_destination, "w", newline="\n") as normal_out, \
            open(tumor_destination, "w", newline="\n") as tumor_out:
        reader = csv.reader(source, delimiter="\t")
        normal_writer = csv.writer(normal_out, delimiter="\t", lineterminator="\n")
        tumor_writer = csv.writer(tumor_out, delimiter="\t", lineterminator="\n")

        try:
            header = next(reader)
        except StopIteration:
            raise RuntimeError(f"{full_mut_tsv} is empty")
        index = {name: i for i, name in enumerate(header)}

        if GERMLINE_COLUMN not in index:
            raise MissingGermlineColumnError(
                f"{full_mut_tsv} has no {GERMLINE_COLUMN} column, which the current "
                f"`msmutect --from_file` requires in the normal histogram.\n{CONVERTER_HINT}")
        for required in ("CHROMOSOME", "START", "END"):
            if required not in index:
                raise RuntimeError(f"{full_mut_tsv} has no {required} column; is it really a "
                                   f"*.full.mut.tsv? (columns: {', '.join(header[:8])}...)")

        plans = {}
        for prefix, writer in ((NORMAL, normal_writer), (TUMOR, tumor_writer)):
            plan = column_plan(index, LOCUS_COLUMNS + histogram_columns(prefix))
            if not any(name.startswith(f"{prefix}MOTIF_REPEATS") for name, _ in plan):
                raise RuntimeError(f"{full_mut_tsv} has no {prefix}MOTIF_REPEATS_* columns, so "
                                   f"no {prefix.rstrip('_').lower()} histogram can be split out")
            writer.writerow([name for name, _ in plan])
            plans[prefix] = (plan, writer)

        normal_plan, _ = plans[NORMAL]
        tumor_plan, _ = plans[TUMOR]
        rows = 0
        for row in reader:
            rows += 1
            normal_writer.writerow([row[i] for _, i in normal_plan])
            tumor_writer.writerow([row[i] for _, i in tumor_plan])
            if progress and rows % PROGRESS_EVERY == 0:
                print(f"    {rows:,} loci...", flush=True)

    if rows == 0:
        raise RuntimeError(f"{full_mut_tsv} has a header but no data rows")
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--full-mut-tsv", required=True, help="*.full.mut.tsv (or .gz) to split")
    parser.add_argument("--output-normal", required=True, help="pseudo-normal histogram to write")
    parser.add_argument("--output-tumor", required=True, help="pseudo-tumor histogram to write")
    args = parser.parse_args()

    try:
        rows = split_full_mut_tsv(args.full_mut_tsv, args.output_normal, args.output_tumor)
    except (MissingGermlineColumnError, RuntimeError, FileNotFoundError) as e:
        sys.exit(f"ERROR: {e}")
    print(f"split {rows:,} loci")
    print(f"  normal -> {args.output_normal}")
    print(f"  tumor  -> {args.output_tumor}")


if __name__ == "__main__":
    main()
