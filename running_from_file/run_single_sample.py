"""Layer 2 of the pyramid: run `msmutect --from_file` on ONE sample.

Two entry points, depending on what you already have:

    run_from_full_mut_tsv(...)   a *.full.mut.tsv -> splits it with layer 1, then runs
    run_from_histograms(...)     already-split normal/tumor histograms -> runs directly

Layer 1 (`split_histograms`) does the splitting; layer 3 (`run_on_directory`) calls this
module once per file. Nothing here knows about directories or sample sets.

The input must already carry KNOWN_GERMLINE_VARIANT -- see split_histograms, and add it
with converting_from_old_format_to_new if it is missing.

Needs MSMuTect's own environment (numpy/scipy/pysam) on PATH, because msmutect.sh runs
whichever `python3` it finds -- e.g. `conda activate genomics`.
"""

import argparse
import csv
import os
import shutil
import subprocess
import sys

from running_from_file import split_histograms

MSMUTECT_CANDIDATES = [
    os.environ.get("MSMUTECT_PATH", ""),
    os.path.expanduser("~/MaruvkaLab/MSMuTect/msmutect.sh"),
]


def find_msmutect(explicit: str = None) -> str:
    if explicit:
        if not os.path.exists(explicit):
            raise FileNotFoundError(f"--msmutect path does not exist: {explicit}")
        return explicit
    for candidate in MSMUTECT_CANDIDATES:
        if candidate and os.path.exists(candidate):
            return candidate
    on_path = shutil.which("msmutect.sh")
    if on_path:
        return on_path
    raise FileNotFoundError(
        "could not find msmutect.sh. Pass --msmutect /path/to/msmutect.sh, or set the "
        "MSMUTECT_PATH environment variable.")


def run_from_histograms(normal_histogram: str, tumor_histogram: str, output_prefix: str,
                        msmutect: str = None, read_level: int = None, vcf: bool = False,
                        quiet: bool = False) -> str:
    """Run `msmutect --from_file` on an already-split pair. Returns the .full.mut.tsv path."""
    msmutect = msmutect or find_msmutect()
    os.makedirs(os.path.dirname(os.path.abspath(output_prefix)), exist_ok=True)

    command = [msmutect, "-N", normal_histogram, "-T", tumor_histogram,
               "-m", "-A", "-O", output_prefix, "--from_file", "-f"]
    if read_level is not None:
        command += ["-r", str(read_level)]
    if vcf:
        command.append("--vcf")

    if not quiet:
        print("  " + " ".join(command), flush=True)
    # captured so a failure inside a worker process still reports usefully
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"msmutect exited {result.returncode} for {normal_histogram}\n"
            f"If that was a ModuleNotFoundError (pysam/numpy/scipy): msmutect.sh runs whichever "
            f"`python3` is on PATH, so activate MSMuTect's environment first "
            f"(e.g. `conda activate genomics`).\n"
            f"command: {' '.join(command)}\n"
            f"stderr:\n{result.stderr[-3000:]}")

    produced = output_prefix + ".full.mut.tsv"
    if not os.path.exists(produced):
        raise RuntimeError(f"msmutect reported success but {produced} was not created")
    return produced


def run_from_full_mut_tsv(full_mut_tsv: str, output_prefix: str, msmutect: str = None,
                          keep_histograms: bool = False, read_level: int = None,
                          vcf: bool = False, quiet: bool = False) -> str:
    """Split `full_mut_tsv` (layer 1) and run msmutect on the result."""
    if os.path.abspath(output_prefix + ".full.mut.tsv") == os.path.abspath(full_mut_tsv):
        raise ValueError(f"output prefix {output_prefix!r} would make msmutect overwrite its own "
                         f"input ({full_mut_tsv}). Choose a different prefix.")
    # resolve msmutect BEFORE the long split, so a bad path fails in a second not an hour
    msmutect = msmutect or find_msmutect()

    normal_histogram = output_prefix + ".pseudo_normal.hist.tsv"
    tumor_histogram = output_prefix + ".pseudo_tumor.hist.tsv"
    rows = split_histograms.split_full_mut_tsv(full_mut_tsv, normal_histogram, tumor_histogram,
                                               progress=not quiet)
    if not quiet:
        print(f"  split {rows:,} loci", flush=True)

    try:
        produced = run_from_histograms(normal_histogram, tumor_histogram, output_prefix,
                                       msmutect=msmutect, read_level=read_level, vcf=vcf,
                                       quiet=quiet)
    except Exception:
        # leave the histograms behind so the run can be retried by hand
        raise
    if not keep_histograms:
        for path in (normal_histogram, tumor_histogram):
            os.remove(path)
    return produced


def call_counts(output_prefix: str) -> dict:
    counts_file = output_prefix + ".full.call_counts.tsv"
    if not os.path.exists(counts_file):
        return {}
    with open(counts_file, newline="") as f:
        return {row["CALL"]: int(row["count"]) for row in csv.DictReader(f, delimiter="\t")}


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    source = parser.add_argument_group("input (give EITHER --full-mut-tsv OR both histograms)")
    source.add_argument("--full-mut-tsv", help="*.full.mut.tsv (or .gz); split, then run")
    source.add_argument("--normal-histogram", help="already-split pseudo-normal histogram")
    source.add_argument("--tumor-histogram", help="already-split pseudo-tumor histogram")
    parser.add_argument("--output-prefix", required=True,
                        help="prefix for everything written, e.g. results/sample_recalled")
    parser.add_argument("--msmutect", default=None,
                        help="path to msmutect.sh (default: $MSMUTECT_PATH, "
                             "~/MaruvkaLab/MSMuTect/msmutect.sh, then PATH)")
    parser.add_argument("--keep-histograms", action="store_true",
                        help="keep the split histograms (several GB genome-wide; deleted on "
                             "success by default). Ignored when passing histograms in.")
    parser.add_argument("--read-level", type=int, default=None,
                        help="msmutect -r: minimum reads to call an allele")
    parser.add_argument("--vcf", action="store_true", help="also ask msmutect for a VCF")
    args = parser.parse_args()

    gave_histograms = bool(args.normal_histogram) and bool(args.tumor_histogram)
    if bool(args.full_mut_tsv) == gave_histograms:
        parser.error("give EITHER --full-mut-tsv OR both --normal-histogram and "
                     "--tumor-histogram, not both and not neither")
    if (bool(args.normal_histogram) != bool(args.tumor_histogram)):
        parser.error("--normal-histogram and --tumor-histogram must be given together")

    try:
        if args.full_mut_tsv:
            produced = run_from_full_mut_tsv(args.full_mut_tsv, args.output_prefix,
                                             msmutect=args.msmutect,
                                             keep_histograms=args.keep_histograms,
                                             read_level=args.read_level, vcf=args.vcf)
        else:
            produced = run_from_histograms(args.normal_histogram, args.tumor_histogram,
                                           args.output_prefix, msmutect=args.msmutect,
                                           read_level=args.read_level, vcf=args.vcf)
    except (split_histograms.MissingGermlineColumnError, RuntimeError,
            FileNotFoundError, ValueError) as e:
        sys.exit(f"ERROR: {e}")

    counts = call_counts(args.output_prefix)
    if counts:
        print("\ncall counts:")
        for call, count in counts.items():
            print(f"  {call:<34}{count:>12,}")
    print(f"\ndone -> {produced}")


if __name__ == "__main__":
    main()
