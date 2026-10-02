#!/usr/bin/env python3
"""Re-call an existing MSMuTect *.full.mut.tsv with the current MSMuTect, in one step.

A *.full.mut.tsv already contains the per-locus normal and tumor histograms, so an old run
can be re-called with the current version without going back to the BAMs. This script does
the whole job:

    1. carves a pseudo-normal and a pseudo-tumor histogram out of the input
    2. adds the KNOWN_GERMLINE_VARIANT column (older MSMuTect predates it, and the current
       --from_file run refuses to start without it in the normal histogram)
    3. runs `msmutect.sh -N ... -T ... -m -A -O <prefix> --from_file`
    4. prints the resulting call counts

The germline annotation is read line-for-line from the loci file rather than looked up,
because a *.full.mut.tsv holds exactly one row per locus in the loci file's order. That
assumption is checked on EVERY row by comparing (CHROMOSOME, START, END); on the first
disagreement the script stops rather than mis-annotate the rest of the genome.

Standalone and stdlib-only -- copy this file anywhere and run it.

Example
-------
    python msmutect_from_file.py \\
        --full-mut-tsv  sample.full.mut.tsv.gz \\
        --loci-file     GRCh38.d1.vd1_1to15_repetitive_loci_sorted_fixed.with_germline_variance \\
        --output-prefix results/sample_recalled

produces results/sample_recalled.full.mut.tsv (plus .mutated_only.mut.tsv and
.call_counts.tsv). Use --carve-only to stop after step 2.

Requires the environment MSMuTect itself needs (numpy, scipy, pysam) -- e.g.
`conda activate genomics`.
"""

import argparse
import csv
import gzip
import os
import shutil
import subprocess
import sys

# ---------------------------------------------------------------- carving

# Locus columns for the carved histograms, in output order. KNOWN_GERMLINE_VARIANT sits
# right after REFERENCE_REPEATS, matching what current MSMuTect writes.
LOCUS_COLUMNS = ["CHROMOSOME", "START", "END", "PATTERN",
                 "REFERENCE_SEQUENCE", "REFERENCE_REPEATS", "KNOWN_GERMLINE_VARIANT"]
GERMLINE_COLUMN = "KNOWN_GERMLINE_VARIANT"

# MSMuTect parses this with str_to_bool, which tests `value.strip() == "True"`, so it has to
# be written as a python bool string -- 1/0 or true/false would silently read as False.
TRUE, FALSE = "True", "False"

HETEROGENEOUS = "HETEROGENEOUS"
VALID_GERMLINE_VALUES = {"HOMOGENEOUS", HETEROGENEOUS}

# 0-indexed fields in the loci file. The germline annotation is its LAST field, read
# positionally off the end so the long sequence fields are never split apart.
LOCI_CHROMOSOME, LOCI_START, LOCI_END = 0, 3, 4

PROGRESS_EVERY = 2_000_000

MSMUTECT_CANDIDATES = [
    os.environ.get("MSMUTECT_PATH", ""),
    os.path.expanduser("~/MaruvkaLab/MSMuTect/msmutect.sh"),
]


def fail(message: str):
    sys.exit(f"ERROR: {message}")


def open_maybe_gzip(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", newline="")
    return open(path, newline="")


def normalize_chromosome(chromosome: str) -> str:
    # loci files use "chr1", mutation files usually use "1"; compare them stripped
    return chromosome[3:] if chromosome.startswith("chr") else chromosome


def histogram_columns(prefix: str):
    return ([f"{prefix}MOTIF_REPEATS_{i}" for i in range(1, 7)] +
            [f"{prefix}SUPPORTING_READS_{i}" for i in range(1, 7)])


def parse_loci_line(line: str, line_number: int):
    """-> ((chromosome, start, end), "HOMOGENEOUS"/"HETEROGENEOUS")"""
    head = line.split("\t", LOCI_END + 1)
    if len(head) <= LOCI_END + 1:
        fail(f"line {line_number} of the loci file has only {len(head)} fields; expected at "
             f"least {LOCI_END + 2} (chromosome, start, end, ... , annotation)")
    annotation = head[-1].rstrip("\n").rsplit("\t", 1)[-1]
    if annotation not in VALID_GERMLINE_VALUES:
        fail(f"line {line_number} of the loci file ends with {annotation!r}, expected one of "
             f"{sorted(VALID_GERMLINE_VALUES)}. Is this a loci file WITH the germline-variance "
             f"column (...with_germline_variance)?")
    return (normalize_chromosome(head[LOCI_CHROMOSOME]), head[LOCI_START], head[LOCI_END]), annotation


def column_plan(header_index: dict, wanted_columns):
    """-> [(output_name, source_index_or_None)]; None means "fill in the germline value"."""
    plan = []
    for name in wanted_columns:
        if name in header_index:
            plan.append((name, header_index[name]))
        elif name == GERMLINE_COLUMN:
            plan.append((name, None))
    return plan


def carve_histograms(full_mut_tsv_path: str, loci_file: str,
                     output_tumor_path: str, output_normal_path: str) -> int:
    """Write both histograms, annotated from the loci file. Returns the row count.

    Exits with an error if the mutation file and the loci file ever disagree on
    (CHROMOSOME, START, END).
    """
    for path, what in ((full_mut_tsv_path, "mutation file"), (loci_file, "loci file")):
        if not os.path.exists(path):
            fail(f"{what} not found: {path}")
    for path in (output_tumor_path, output_normal_path):
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)

    with open_maybe_gzip(full_mut_tsv_path) as source, \
            open(loci_file) as loci, \
            open(output_normal_path, "w", newline="\n") as normal_out, \
            open(output_tumor_path, "w", newline="\n") as tumor_out:

        reader = csv.reader(source, delimiter="\t")
        normal_writer = csv.writer(normal_out, delimiter="\t", lineterminator="\n")
        tumor_writer = csv.writer(tumor_out, delimiter="\t", lineterminator="\n")

        try:
            header = next(reader)
        except StopIteration:
            fail(f"{full_mut_tsv_path} is empty")
        index = {name: i for i, name in enumerate(header)}
        for required in ("CHROMOSOME", "START", "END"):
            if required not in index:
                fail(f"{full_mut_tsv_path} has no {required} column; is it really a "
                     f"*.full.mut.tsv? (columns: {', '.join(header[:8])}...)")

        normal_plan = column_plan(index, LOCUS_COLUMNS + histogram_columns("NORMAL_"))
        tumor_plan = column_plan(index, LOCUS_COLUMNS + histogram_columns("TUMOR_"))
        for plan, prefix, path in ((normal_plan, "NORMAL_", output_normal_path),
                                   (tumor_plan, "TUMOR_", output_tumor_path)):
            if not any(name.startswith(f"{prefix}MOTIF_REPEATS") for name, _ in plan):
                fail(f"{full_mut_tsv_path} has no {prefix}MOTIF_REPEATS_* columns, so no "
                     f"{prefix.rstrip('_').lower()} histogram can be carved for {path}")

        if GERMLINE_COLUMN in index:
            print(f"note: input already has {GERMLINE_COLUMN}; it is being REPLACED with the "
                  f"value from {os.path.basename(loci_file)}")
            normal_plan = [(n, None if n == GERMLINE_COLUMN else i) for n, i in normal_plan]
            tumor_plan = [(n, None if n == GERMLINE_COLUMN else i) for n, i in tumor_plan]

        normal_writer.writerow([name for name, _ in normal_plan])
        tumor_writer.writerow([name for name, _ in tumor_plan])

        i_chromosome, i_start, i_end = (index[c] for c in ("CHROMOSOME", "START", "END"))
        rows = germline_count = 0
        for row in reader:
            rows += 1
            loci_line = loci.readline()
            if not loci_line:
                fail(f"the loci file ran out after {rows - 1} loci, but {full_mut_tsv_path} "
                     f"still has rows (row {rows}). The two files do not describe the same "
                     f"loci set.")
            loci_locus, annotation = parse_loci_line(loci_line, rows)
            mut_locus = (normalize_chromosome(row[i_chromosome]), row[i_start], row[i_end])
            if mut_locus != loci_locus:
                fail(f"locus mismatch at row {rows}: {os.path.basename(full_mut_tsv_path)} has "
                     f"{mut_locus}, {os.path.basename(loci_file)} has {loci_locus}.\n"
                     f"       The germline annotation is taken line-for-line, so the files "
                     f"must list the same loci in the same order. Nothing further was written "
                     f"-- check that this is the loci file the run used.")

            germline = TRUE if annotation == HETEROGENEOUS else FALSE
            germline_count += annotation == HETEROGENEOUS
            normal_writer.writerow([germline if i is None else row[i] for _, i in normal_plan])
            tumor_writer.writerow([germline if i is None else row[i] for _, i in tumor_plan])
            if rows % PROGRESS_EVERY == 0:
                print(f"  {rows:,} loci...", flush=True)

        leftover_loci = sum(1 for _ in loci)

    if rows == 0:
        fail(f"{full_mut_tsv_path} has a header but no data rows")

    print(f"carved {rows:,} loci "
          f"({germline_count:,} / {100 * germline_count / rows:.1f}% {GERMLINE_COLUMN}=True)")
    print(f"  normal histogram -> {output_normal_path}")
    print(f"  tumor  histogram -> {output_tumor_path}")
    if leftover_loci:
        # every row written was verified aligned, so the output is usable; but the input
        # covered only part of the loci file, which is worth saying out loud
        print(f"WARNING: the loci file has {leftover_loci:,} more loci than the mutation file. "
              f"All {rows:,} carved loci matched, but this input covers only part of it.")
    return rows


# ---------------------------------------------------------------- msmutect


def find_msmutect(explicit: str = None) -> str:
    if explicit:
        if not os.path.exists(explicit):
            fail(f"--msmutect path does not exist: {explicit}")
        return explicit
    for candidate in MSMUTECT_CANDIDATES:
        if candidate and os.path.exists(candidate):
            return candidate
    on_path = shutil.which("msmutect.sh")
    if on_path:
        return on_path
    fail("could not find msmutect.sh. Pass --msmutect /path/to/msmutect.sh, or set "
         "the MSMUTECT_PATH environment variable.")


def run_msmutect(msmutect: str, normal_hist: str, tumor_hist: str, output_prefix: str,
                 read_level: int = None, vcf: bool = False) -> str:
    command = [msmutect, "-N", normal_hist, "-T", tumor_hist,
               "-m", "-A", "-O", output_prefix, "--from_file", "-f"]
    if read_level is not None:
        command += ["-r", str(read_level)]
    if vcf:
        command.append("--vcf")

    print("\nrunning: " + " ".join(command), flush=True)
    result = subprocess.run(command)
    if result.returncode != 0:
        fail(f"msmutect exited {result.returncode}.\n"
             f"       If that was a ModuleNotFoundError (pysam/numpy/scipy): msmutect.sh runs\n"
             f"       whichever `python3` is on PATH, so activate MSMuTect's environment first\n"
             f"       (e.g. `conda activate genomics`) and re-run.\n"
             f"       The histograms were kept, so you can also just run it by hand:\n"
             f"       {' '.join(command)}")

    produced = output_prefix + ".full.mut.tsv"
    if not os.path.exists(produced):
        fail(f"msmutect reported success but {produced} was not created")
    return produced


def print_call_counts(output_prefix: str):
    counts_file = output_prefix + ".full.call_counts.tsv"
    if not os.path.exists(counts_file):
        return
    print("\ncall counts:")
    with open(counts_file, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            print(f"  {row['CALL']:<34}{int(row['count']):>12,}")


# ---------------------------------------------------------------- cli


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--full-mut-tsv", required=True,
                        help="existing *.full.mut.tsv (or .gz) to re-call")
    parser.add_argument("--loci-file", required=True,
                        help="loci file whose last column is the HOMOGENEOUS/HETEROGENEOUS "
                             "germline-variance annotation")
    parser.add_argument("--output-prefix", required=True,
                        help="prefix for everything written, e.g. results/sample_recalled")
    parser.add_argument("--msmutect", default=None,
                        help="path to msmutect.sh (default: $MSMUTECT_PATH, "
                             "~/MaruvkaLab/MSMuTect/msmutect.sh, then PATH)")
    parser.add_argument("--carve-only", action="store_true",
                        help="write the histograms and stop, without running msmutect")
    parser.add_argument("--keep-histograms", action="store_true",
                        help="keep the carved histograms (they are several GB for a whole "
                             "genome; deleted after a successful run by default)")
    parser.add_argument("--read-level", type=int, default=None,
                        help="msmutect -r: minimum reads to call an allele")
    parser.add_argument("--vcf", action="store_true", help="also ask msmutect for a VCF")
    args = parser.parse_args()

    normal_hist = args.output_prefix + ".pseudo_normal.hist.tsv"
    tumor_hist = args.output_prefix + ".pseudo_tumor.hist.tsv"
    would_produce = os.path.abspath(args.output_prefix + ".full.mut.tsv")
    if would_produce == os.path.abspath(args.full_mut_tsv):
        fail(f"--output-prefix would make msmutect overwrite the input file "
             f"({args.full_mut_tsv}). Choose a different prefix.")

    # resolve msmutect BEFORE the long carve, so a bad path fails in a second not an hour
    msmutect = None if args.carve_only else find_msmutect(args.msmutect)

    carve_histograms(args.full_mut_tsv, args.loci_file, tumor_hist, normal_hist)
    if args.carve_only:
        print("\n--carve-only: stopping before msmutect. Run it with:")
        print(f"  msmutect.sh -N {normal_hist} -T {tumor_hist} "
              f"-m -A -O {args.output_prefix} --from_file -f")
        return

    produced = run_msmutect(msmutect, normal_hist, tumor_hist, args.output_prefix,
                            read_level=args.read_level, vcf=args.vcf)
    print_call_counts(args.output_prefix)

    if not args.keep_histograms:
        for path in (normal_hist, tumor_hist):
            os.remove(path)
        print("\n(carved histograms deleted; --keep-histograms to retain them)")

    print(f"\ndone -> {produced}")
    for suffix in (".full.mutated_only.mut.tsv", ".full.call_counts.tsv"):
        extra = args.output_prefix + suffix
        if os.path.exists(extra):
            print(f"        {extra}")


if __name__ == "__main__":
    main()
