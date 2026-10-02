"""Re-run MSMuTect from file over the GIAB sample set.

A thin caller of layer 3 (`run_on_directory`) with the GIAB paths filled in. For same-sample
GIAB pairs the resulting mutation calls are false positives, so this measures the FPR.

This used to carve the histograms itself with HARDCODED column positions:

    cut -f1-19       IN > <name>.pseudo_normal.hist.tsv
    cut -f1-6,29-41  IN > <name>.pseudo_tumor.hist.tsv

which silently produces the wrong histograms the moment a column is added or reordered --
and it could not add KNOWN_GERMLINE_VARIANT at all, so it cannot feed current MSMuTect. The
pyramid selects columns by NAME instead, so it survives schema changes.

Inputs must already carry KNOWN_GERMLINE_VARIANT; add it with
converting_from_old_format_to_new/convert_old_msmutect_file_to_new.py.
"""

import argparse

from running_from_file import run_on_directory

INPUT_DIR = "/data/gib_files_filtered"
OUTPUT_DIR = "/results/filtered_through_msmutect_directly"
PATTERN = "*.filtered.full.mut.tsv"


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input-dir", default=INPUT_DIR)
    parser.add_argument("--output-dir", default=OUTPUT_DIR)
    parser.add_argument("--pattern", nargs="+", default=[PATTERN])
    parser.add_argument("--workers", type=int, default=run_on_directory.DEFAULT_WORKERS)
    parser.add_argument("--force", action="store_true", help="re-run even if the output exists")
    args = parser.parse_args()

    results = run_on_directory.run_directory(args.input_dir, args.output_dir,
                                             patterns=tuple(args.pattern),
                                             workers=args.workers, force=args.force)
    if any(produced is None for produced, _ in results.values()):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
