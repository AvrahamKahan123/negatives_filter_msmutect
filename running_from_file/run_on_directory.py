"""Layer 3 of the pyramid: run `msmutect --from_file` over a whole directory.

Globs a directory for *.full.mut.tsv files and hands each one to layer 2
(`run_single_sample.run_from_full_mut_tsv`), which in turn uses layer 1 to split it.
This module adds only the fan-out: discovery, parallelism, skip-existing and reporting.

Each --from_file run is single-core by design, so the way to use the machine is to run
several samples at once rather than to try to speed up one.

    python -m running_from_file.run_on_directory \\
        --input-dir  data/filtered_full_old_wgermline \\
        --output-dir results/recalled

Inputs must already carry KNOWN_GERMLINE_VARIANT (see converting_from_old_format_to_new);
a file missing it is reported and skipped rather than killing the whole run.
"""

import argparse
import glob
import multiprocessing as mp
import os
import time

from running_from_file import run_single_sample, split_histograms

DEFAULT_PATTERNS = ("*.full.mut.tsv", "*.full.mut.tsv.gz", "*.filt.mut.tsv", "*.filt.mut.tsv.gz")
DEFAULT_WORKERS = 6


def input_files(input_dir: str, patterns=DEFAULT_PATTERNS):
    found = []
    for pattern in patterns:
        found.extend(glob.glob(os.path.join(input_dir, pattern)))
    # a sample present as both .gz and plain would otherwise be run twice
    by_sample = {}
    for path in sorted(found):
        by_sample.setdefault(split_histograms.sample_name(path), path)
    return sorted(by_sample.items())


def _run_one(job) -> tuple:
    """(sample, input) -> (sample, output_or_None, error_or_None). Never raises, so one bad
    file cannot take down the pool."""
    sample, source, output_dir, msmutect, read_level, vcf, force = job
    output_prefix = os.path.join(output_dir, sample)
    produced = output_prefix + ".full.mut.tsv"
    if os.path.exists(produced) and not force:
        return sample, produced, "skipped (already exists)"
    try:
        out = run_single_sample.run_from_full_mut_tsv(
            source, output_prefix, msmutect=msmutect, read_level=read_level, vcf=vcf, quiet=True)
        return sample, out, None
    except split_histograms.MissingGermlineColumnError as e:
        return sample, None, f"needs conversion: {e}"
    except Exception as e:  # noqa: BLE001 - a worker must report, not crash the pool
        return sample, None, str(e)


def run_directory(input_dir: str, output_dir: str, patterns=DEFAULT_PATTERNS,
                  workers: int = DEFAULT_WORKERS, msmutect: str = None,
                  read_level: int = None, vcf: bool = False, force: bool = False) -> dict:
    files = input_files(input_dir, patterns)
    if not files:
        raise SystemExit(f"no files matching {', '.join(patterns)} in {input_dir}")
    os.makedirs(output_dir, exist_ok=True)
    msmutect = run_single_sample.find_msmutect(msmutect)  # fail fast, before any work

    jobs = [(sample, source, output_dir, msmutect, read_level, vcf, force)
            for sample, source in files]
    print(f"{len(jobs)} file(s) from {input_dir}, {workers} at a time -> {output_dir}")

    started = time.time()
    results = {}
    with mp.Pool(processes=workers) as pool:
        for i, (sample, produced, note) in enumerate(pool.imap_unordered(_run_one, jobs), 1):
            elapsed = (time.time() - started) / 60
            if produced and not note:
                status = f"ok -> {os.path.basename(produced)}"
            elif produced:
                status = note
            else:
                status = f"FAILED: {note}"
            print(f"[{i}/{len(jobs)}] {sample:<28} {status}   [{elapsed:.1f} min]", flush=True)
            results[sample] = (produced, note)

    failed = [s for s, (produced, _) in results.items() if produced is None]
    print(f"\ndone in {(time.time() - started) / 60:.1f} min: "
          f"{len(results) - len(failed)}/{len(results)} succeeded -> {output_dir}")
    if failed:
        print(f"FAILED ({len(failed)}): {', '.join(failed)}")
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--pattern", nargs="+", default=list(DEFAULT_PATTERNS),
                        help=f"glob(s) to match (default: {' '.join(DEFAULT_PATTERNS)})")
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS,
                        help=f"concurrent msmutect runs (default: {DEFAULT_WORKERS})")
    parser.add_argument("--msmutect", default=None, help="path to msmutect.sh")
    parser.add_argument("--read-level", type=int, default=None, help="msmutect -r")
    parser.add_argument("--vcf", action="store_true", help="also ask msmutect for a VCF")
    parser.add_argument("--force", action="store_true", help="re-run even if the output exists")
    args = parser.parse_args()

    results = run_directory(args.input_dir, args.output_dir, patterns=tuple(args.pattern),
                            workers=args.workers, msmutect=args.msmutect,
                            read_level=args.read_level, vcf=args.vcf, force=args.force)
    if any(produced is None for produced, _ in results.values()):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
