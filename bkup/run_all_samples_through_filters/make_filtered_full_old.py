import argparse
import glob
import os

import config
import old_mutation_calls

# Writes a mutation-only copy of every *.full.mut.tsv[.gz] in data/full/ to
# data/filtered_full_old/, one output file per input:
#
#     data/full/TCGA-A6-5661.full.mut.tsv.gz
#         -> data/filtered_full_old/TCGA-A6-5661.filtered.full.mut.tsv
#
# "Mutation lines" means the rows the old run itself called CALL == "M". The header is
# kept, so each output is a valid *.mut.tsv on its own. Everything else -- NM, AN, INS, RR,
# FFT, TMA, LOH -- is dropped, which is ~27M of the ~27.7M rows.
#
# Note these come out roughly 300 MB each uncompressed (~1.5 GB for all 5); pass --gzip to
# write *.filtered.full.mut.tsv.gz instead.

FULL_SUFFIXES = (".full.mut.tsv.gz", ".full.mut.tsv")


def sample_name(path: str) -> str:
    name = os.path.basename(path)
    for suffix in FULL_SUFFIXES:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return os.path.splitext(name)[0]


def input_files(samples=None):
    found = []
    for suffix in FULL_SUFFIXES:
        found.extend(glob.glob(os.path.join(config.FULL_DIR, f"*{suffix}")))
    # a sample present as both .gz and plain would otherwise be processed twice
    by_sample = {}
    for path in sorted(found):
        by_sample.setdefault(sample_name(path), path)
    if samples:
        missing = [s for s in samples if s not in by_sample]
        if missing:
            raise SystemExit(f"no full file in {config.FULL_DIR} for: {', '.join(missing)}")
        by_sample = {s: by_sample[s] for s in samples}
    return sorted(by_sample.items())


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--samples", nargs="+", default=None,
                         help="only these samples (default: every full file in data/full/)")
    parser.add_argument("--force", action="store_true", help="rewrite even if the output exists")
    parser.add_argument("--gzip", action="store_true", help="write the outputs gzipped")
    args = parser.parse_args()

    files = input_files(args.samples)
    if not files:
        raise SystemExit(f"no *.full.mut.tsv[.gz] found in {config.FULL_DIR} -- "
                         f"run download_all.py first")

    os.makedirs(config.FILTERED_FULL_OLD_DIR, exist_ok=True)
    config.ensure_free_disk_gb()
    print(f"{len(files)} file(s) -> {config.FILTERED_FULL_OLD_DIR}")

    total_kept = 0
    for i, (sample, source) in enumerate(files, 1):
        destination = config.filtered_full_old_path(sample)
        if args.gzip:
            destination += ".gz"
        if os.path.exists(destination) and not args.force:
            print(f"[{i}/{len(files)}] {sample}: already exists, skipping (--force to rewrite)")
            continue

        print(f"[{i}/{len(files)}] {sample}: {os.path.basename(source)}")
        kept, total = old_mutation_calls.write_mutation_calls(source, destination)
        total_kept += kept
        size_gb = os.path.getsize(destination) / 1024 ** 3
        print(f"    {kept:,} mutations kept of {total:,} loci "
              f"({100 * kept / max(total, 1):.2f}%) -> {destination} [{size_gb:.2f} GB]")

    print(f"\ndone: {total_kept:,} mutation rows written to {config.FILTERED_FULL_OLD_DIR}")


if __name__ == "__main__":
    main()
