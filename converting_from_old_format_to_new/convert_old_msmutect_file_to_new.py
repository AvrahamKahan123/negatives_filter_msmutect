"""Convert old-format msmutect output to the new format.

The two "Noisy Locus" columns are dropped and a KNOWN_GERMLINE_VARIANT column is inserted
directly after REFERENCE_REPEATS, populated from the noisy-locus / known-germline-variant
database. All other columns are preserved.

    # one file
    python -m converting_from_old_format_to_new.convert_old_msmutect_file_to_new \\
        --input old.full.mut.tsv --output new.full.mut.tsv

    # a whole directory
    python -m converting_from_old_format_to_new.convert_old_msmutect_file_to_new \\
        --input-dir data/old_format --output-dir data/old_format_wgermline

The result is what `running_from_file` expects: its layer 1 refuses input without
KNOWN_GERMLINE_VARIANT, because adding it is this module's job.
"""

import argparse
import glob
import gzip
import os
import re
import sys
from typing import List, Set

from converting_from_old_format_to_new.NoisyLocusDB import NoisyLocusDB

# Column that existed (twice: once for normal, once for tumor) in the old msmutect
# output and is dropped in the new format.
NOISY_LOCUS_COLUMN = "Noisy Locus"
# New column added by this conversion, populated from the germline/noisy locus DB.
KNOWN_GERMLINE_VARIANT_COLUMN = "KNOWN_GERMLINE_VARIANT"
# In the new format KNOWN_GERMLINE_VARIANT sits immediately after this column.
INSERT_AFTER_COLUMN = "REFERENCE_REPEATS"

# Locus-identifying columns are at the front of the row and are unaffected by the
# column additions/removals, so their indices are stable.
CHROMOSOME_COL = 0
START_COL = 1
END_COL = 2
PATTERN_COL = 3

DEFAULT_PATTERN = "*.tsv"


class ConversionError(RuntimeError):
    """Input is not an old-format msmutect file this module can convert."""


def _open_read(path: str):
    # the CLI accepts arbitrary paths, and the lab's files are routinely gzipped
    if path.endswith(".gz"):
        return gzip.open(path, "rt", newline="")
    return open(path, "r", newline="")


def _open_write(path: str, gzipped: bool):
    # `gzipped` comes from the real destination name, not `path`, because output is written
    # to a .partial temporary whose own extension says nothing about the encoding
    if gzipped:
        return gzip.open(path, "wt", newline="")
    return open(path, "w", newline="")


def _is_noisy_locus_column(name: str) -> bool:
    """True for the noisy-locus columns, including pandas-deduplicated names.

    A file may hold the column twice as "Noisy Locus" and "Noisy Locus.1"
    (pandas renames duplicate headers), so strip any trailing ".<n>" suffix
    before comparing.
    """
    base_name = re.sub(r"\.\d+$", "", name)
    return base_name == NOISY_LOCUS_COLUMN


def _transform_row(fields: List[str], inserted_value: str, drop_indices: Set[int], insert_after: int) -> List[str]:
    """Drop the noisy-locus columns and insert `inserted_value` right after `insert_after`."""
    transformed = []
    for i, value in enumerate(fields):
        if i in drop_indices:
            continue
        transformed.append(value)
        if i == insert_after:
            transformed.append(inserted_value)
    return transformed


def convert_old_msmutect_file_to_new(input_fp: str, output_fp: str,
                                     noisy_db: NoisyLocusDB = None) -> int:
    """Convert an old-format msmutect output file to the new format. Returns the row count.

    The two "Noisy Locus" columns are removed, and a KNOWN_GERMLINE_VARIANT column
    is inserted directly after REFERENCE_REPEATS, populated by querying the noisy
    locus / known-germline-variant database. Any other columns are preserved.

    `noisy_db` may be supplied so a single database (which is expensive to load) can
    be reused across many files; when omitted one is loaded for this file.

    Input and output may be plain or gzipped (.gz). Raises ConversionError if the input
    is not convertible -- in particular if it has ALREADY been converted, which would
    otherwise silently produce a second KNOWN_GERMLINE_VARIANT column.
    """
    if noisy_db is None:
        noisy_db = NoisyLocusDB()

    os.makedirs(os.path.dirname(os.path.abspath(output_fp)), exist_ok=True)
    # write to .partial and rename only on success, so an interrupted run never leaves a
    # truncated file that convert_directory's skip-existing check would happily reuse
    partial_fp = output_fp + ".partial"
    rows = 0

    with _open_read(input_fp) as in_file:
        # validate the header BEFORE creating any output, so a rejected input leaves nothing
        header_line = in_file.readline()
        if not header_line:
            raise ConversionError(f"{input_fp} is empty")
        header = header_line.rstrip("\n").split("\t")

        if KNOWN_GERMLINE_VARIANT_COLUMN in header:
            raise ConversionError(
                f"{input_fp} already has a {KNOWN_GERMLINE_VARIANT_COLUMN} column -- it looks "
                f"like it has already been converted. Converting again would give it two.")
        if INSERT_AFTER_COLUMN not in header:
            raise ConversionError(
                f"{input_fp} has no {INSERT_AFTER_COLUMN} column, so there is nowhere to insert "
                f"{KNOWN_GERMLINE_VARIANT_COLUMN}. Is it really a msmutect output file? "
                f"(columns: {', '.join(header[:8])}...)")

        drop_indices = {i for i, name in enumerate(header) if _is_noisy_locus_column(name)}
        insert_after = header.index(INSERT_AFTER_COLUMN)

        try:
            with _open_write(partial_fp, output_fp.endswith(".gz")) as out_file:
                new_header = _transform_row(header, KNOWN_GERMLINE_VARIANT_COLUMN,
                                            drop_indices, insert_after)
                out_file.write("\t".join(new_header) + "\n")

                for line in in_file:
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) <= PATTERN_COL:
                        raise ConversionError(
                            f"{input_fp} row {rows + 1} has only {len(fields)} column(s); the "
                            f"first {PATTERN_COL + 1} must be CHROMOSOME, START, END, PATTERN. "
                            f"Is the file truncated or not tab-separated?")
                    is_known_germline_variant = noisy_db.locus_present_in_db(
                        fields[CHROMOSOME_COL], fields[START_COL], fields[END_COL], fields[PATTERN_COL]
                    )
                    new_row = _transform_row(fields, str(is_known_germline_variant),
                                             drop_indices, insert_after)
                    out_file.write("\t".join(new_row) + "\n")
                    rows += 1
        except BaseException:
            # never leave a half-written .partial lying next to the real outputs
            if os.path.exists(partial_fp):
                os.remove(partial_fp)
            raise

    os.replace(partial_fp, output_fp)
    return rows


def convert_directory(input_dir: str, output_dir: str, pattern: str = DEFAULT_PATTERN,
                      force: bool = False, noisy_db: NoisyLocusDB = None) -> dict:
    """Convert every file matching `pattern` in `input_dir`, writing same-named results to
    `output_dir`. Returns {filename: rows or None-on-failure}.

    One database is loaded and reused across all files. A file that cannot be converted is
    reported and skipped rather than aborting the whole directory.
    """
    input_files = sorted(glob.glob(os.path.join(input_dir, pattern)))
    if not input_files:
        raise SystemExit(f"no files matching {pattern!r} in {input_dir}")
    os.makedirs(output_dir, exist_ok=True)

    if noisy_db is None:
        print("loading NoisyLocusDB...")
        noisy_db = NoisyLocusDB()
        print(f"loaded {len(noisy_db.db):,} loci")

    results = {}
    for i, input_fp in enumerate(input_files, start=1):
        name = os.path.basename(input_fp)
        output_fp = os.path.join(output_dir, name)
        if os.path.exists(output_fp) and not force:
            print(f"[{i}/{len(input_files)}] {name}: already converted, skipping (--force to redo)")
            results[name] = None
            continue
        print(f"[{i}/{len(input_files)}] {name}", flush=True)
        try:
            rows = convert_old_msmutect_file_to_new(input_fp, output_fp, noisy_db=noisy_db)
            print(f"    {rows:,} rows -> {output_fp}")
            results[name] = rows
        except ConversionError as e:
            print(f"    SKIPPED: {e}")
            results[name] = None
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    single = parser.add_argument_group("one file")
    single.add_argument("--input", help="old-format file to convert (plain or .gz)")
    single.add_argument("--output", help="where to write the converted file (.gz to compress)")
    whole_dir = parser.add_argument_group("a whole directory")
    whole_dir.add_argument("--input-dir", help="directory of old-format files")
    whole_dir.add_argument("--output-dir", help="directory to write converted files to")
    whole_dir.add_argument("--pattern", default=DEFAULT_PATTERN,
                           help=f"glob to match inside --input-dir (default: {DEFAULT_PATTERN})")
    whole_dir.add_argument("--force", action="store_true",
                           help="re-convert even if the output already exists")
    args = parser.parse_args()

    single_mode = bool(args.input) or bool(args.output)
    directory_mode = bool(args.input_dir) or bool(args.output_dir)
    if single_mode == directory_mode:
        parser.error("give EITHER --input/--output OR --input-dir/--output-dir")
    if single_mode and not (args.input and args.output):
        parser.error("--input and --output must be given together")
    if directory_mode and not (args.input_dir and args.output_dir):
        parser.error("--input-dir and --output-dir must be given together")

    if single_mode:
        try:
            rows = convert_old_msmutect_file_to_new(args.input, args.output)
        except (ConversionError, FileNotFoundError) as e:
            sys.exit(f"ERROR: {e}")
        print(f"converted {rows:,} rows -> {args.output}")
        return

    results = convert_directory(args.input_dir, args.output_dir, pattern=args.pattern,
                                force=args.force)
    converted = sum(1 for rows in results.values() if rows is not None)
    print(f"\ndone: {converted}/{len(results)} file(s) converted -> {args.output_dir}")


if __name__ == '__main__':
    main()

