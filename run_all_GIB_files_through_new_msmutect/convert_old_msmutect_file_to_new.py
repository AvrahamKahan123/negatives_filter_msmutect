import os
import re
import glob
from typing import List, Set

from results_postprocessing.NoisyLocusDB import NoisyLocusDB
from results_postprocessing.cancer_data import data_directory

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


def convert_old_msmutect_file_to_new(input_fp: str, output_fp: str, noisy_db: NoisyLocusDB = None):
    """Convert an old-format msmutect output file to the new format.

    The two "Noisy Locus" columns are removed, and a KNOWN_GERMLINE_VARIANT column
    is inserted directly after REFERENCE_REPEATS, populated by querying the noisy
    locus / known-germline-variant database. Any other columns are preserved.

    `noisy_db` may be supplied so a single database (which is expensive to load) can
    be reused across many files; when omitted one is loaded for this file.
    """
    if noisy_db is None:
        noisy_db = NoisyLocusDB()
    with open(input_fp, 'r') as in_file, open(output_fp, 'w') as out_file:
        header = in_file.readline().rstrip("\n").split("\t")
        drop_indices = {i for i, name in enumerate(header) if _is_noisy_locus_column(name)}
        insert_after = header.index(INSERT_AFTER_COLUMN)

        new_header = _transform_row(header, KNOWN_GERMLINE_VARIANT_COLUMN, drop_indices, insert_after)
        out_file.write("\t".join(new_header) + "\n")

        for line in in_file:
            fields = line.rstrip("\n").split("\t")
            is_known_germline_variant = noisy_db.locus_present_in_db(
                fields[CHROMOSOME_COL], fields[START_COL], fields[END_COL], fields[PATTERN_COL]
            )
            new_row = _transform_row(fields, str(is_known_germline_variant), drop_indices, insert_after)
            out_file.write("\t".join(new_row) + "\n")


def convert_directory(input_dir: str, output_dir: str):
    """Convert every .tsv in `input_dir` to the new format, writing each result to
    a file of the same name in `output_dir`."""
    os.makedirs(output_dir, exist_ok=True)
    noisy_db = NoisyLocusDB()
    input_files = sorted(glob.glob(os.path.join(input_dir, "*.tsv")))
    for i, input_fp in enumerate(input_files, start=1):
        output_fp = os.path.join(output_dir, os.path.basename(input_fp))
        print(f"[{i}/{len(input_files)}] {os.path.basename(input_fp)}")
        convert_old_msmutect_file_to_new(input_fp, output_fp, noisy_db=noisy_db)


if __name__ == '__main__':
    convert_directory(os.path.join(data_directory(), "gib_files_filtered"),
                      os.path.join(data_directory(), "gib_files_filtered_new_format"))

