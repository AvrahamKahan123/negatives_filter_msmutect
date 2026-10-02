"""Step 2: add the KNOWN_GERMLINE_VARIANT column.

One condor job per sample:

    <work_root>/filtered/TCGA-W5-AA2K.full.mut.tsv.gz
        -> <work_root>/filtered_wgermline/TCGA-W5-AA2K.full.mut.tsv.gz

All the work is done by converting_from_old_format_to_new: the two "Noisy Locus" columns
are dropped and KNOWN_GERMLINE_VARIANT is inserted after REFERENCE_REPEATS, looked up per
locus in NoisyLocusDB. This module only wires the cluster paths to it.

MEMORY: NoisyLocusDB holds ~3.5M loci in a python set and is loaded once per job, so these
jobs need roughly 2-3 GB -- see step2_germline.sub. That is the reason this step is far
hungrier than the other two.

    python -m run_msmutect_over_whole_tcga_again.add_germline --sample TCGA-W5-AA2K
"""

import argparse
import os

from run_msmutect_over_whole_tcga_again import config  # also puts repo_root on sys.path

from converting_from_old_format_to_new.convert_old_msmutect_file_to_new import (
    ConversionError, convert_old_msmutect_file_to_new)


def add_germline(sample: str, force: bool = False, quiet: bool = False) -> str:
    source = config.filtered_path(sample)
    destination = config.filtered_wgermline_path(sample)
    if not os.path.exists(source):
        raise SystemExit(f"ERROR: {source} not found -- has step 1 run for this sample?")
    if os.path.exists(destination) and not force:
        if not quiet:
            print(f"[{sample}] already annotated, skipping (--force to redo): {destination}")
        return destination

    os.makedirs(os.path.dirname(destination), exist_ok=True)
    try:
        rows = convert_old_msmutect_file_to_new(source, destination)
    except ConversionError as e:
        raise SystemExit(f"ERROR: {e}")
    if not quiet:
        print(f"[{sample}] {rows:,} rows annotated -> {destination}")

    if config.DELETE_FILTERED_AFTER_GERMLINE:
        os.remove(source)
        if not quiet:
            print(f"[{sample}] removed {source} (delete_filtered_after_germline)")
    return destination


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sample", required=True, help="e.g. TCGA-W5-AA2K")
    parser.add_argument("--force", action="store_true", help="re-annotate even if the output exists")
    args = parser.parse_args()
    add_germline(args.sample, force=args.force)


if __name__ == "__main__":
    main()
