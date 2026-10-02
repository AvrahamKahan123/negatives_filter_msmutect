import argparse
import csv
import gzip
import os

import config
import download_samples
from running_from_file import run_single_sample
from running_from_file.split_histograms import histogram_columns  # noqa: F401 (re-exported)

# Step 3 of the Negatives pipeline (Negatives/run_all_samples_through_newest_filter/Plan.md):
# run the newest MSMuTect on each sample's *.full.mut.tsv.gz via `--from_file`.
#
# This is an ADAPTER, not part of the pyramid. The pyramid's layer 1 (split_histograms)
# deliberately refuses input without KNOWN_GERMLINE_VARIANT, because adding that column is
# converting_from_old_format_to_new's job. But this pipeline's inputs are the raw downloaded
# *.full.mut.tsv.gz, which predate the column -- and converting them first would mean an
# extra full read+write of a ~5 GB file per sample. So this module keeps its own carve,
# which adds the column in the SAME pass that splits the histograms, and then hands off to
# layer 2 (run_single_sample.run_from_histograms) for the actual msmutect invocation.
#
# If you are not on that pipeline, use the pyramid instead:
#     convert_old_msmutect_file_to_new.py  ->  run_single_sample / run_on_directory
#
# Histograms are carved by streaming straight out of the .gz (never writing an uncompressed
# copy of the full file to disk), to keep peak disk usage low.
#
# KNOWN_GERMLINE_VARIANT: these files come from an older MSMuTect that predates the column,
# but the current --from_file run REFUSES to start without it in the normal histogram
# (src/Entry/PairFileBatches.py::run_from_file), so we have to add it ourselves.
#
# A *.full.mut.tsv holds one row per locus in exactly the order of the loci file, so rather
# than looking each locus up, we read config.GERMLINE_LOCI_FILE alongside the mutation file
# and take its HOMOGENEOUS/HETEROGENEOUS column line-for-line. Verified exact on a full run:
# 27,705,407 rows, zero locus mismatches, and neither file had a single leftover row.
# Every row still re-checks that the two files agree on the locus, so if a sample ever fails
# to line up the run aborts loudly instead of silently mis-annotating the whole genome.
#
# The value is written as the literal "True"/"False" because MSMuTect parses it with
# str_to_bool (src/Entry/ResultsReaders.py), which tests `value.strip() == "True"`.

GERMLINE_COLUMN = "KNOWN_GERMLINE_VARIANT"
HETEROGENEOUS = "HETEROGENEOUS"
VALID_LOCI_GERMLINE_VALUES = {"HOMOGENEOUS", HETEROGENEOUS}


def parse_loci_line(line: str):
    # returns ((chromosome, start, end, pattern), "HOMOGENEOUS"/"HETEROGENEOUS").
    # splits only up to the pattern field and grabs the annotation off the end, so the long
    # reference/flanking sequence fields are never split into separate strings.
    head = line.split("\t", config.LOCI_FILE_PATTERN + 1)
    annotation = head[-1].rstrip("\n").rsplit("\t", 1)[-1]
    return (head[config.LOCI_FILE_CHROMOSOME][3:],   # strip the "chr" prefix
            head[config.LOCI_FILE_START],
            head[config.LOCI_FILE_END],
            head[config.LOCI_FILE_PATTERN]), annotation


def column_plan(header_index: dict, wanted_columns):
    # (name, source_index) per output column, in wanted_columns order; source_index is None
    # for KNOWN_GERMLINE_VARIANT when the old file lacks it, meaning "compute this one".
    # Columns absent for any other reason are simply skipped.
    plan = []
    for name in wanted_columns:
        if name in header_index:
            plan.append((name, header_index[name]))
        elif name == GERMLINE_COLUMN:
            plan.append((name, None))
    return plan


def carve_histograms(source_gz: str, normal_dest: str, tumor_dest: str):
    # Select the wanted columns BY NAME from the source header and write the pseudo-normal
    # and pseudo-tumor histograms in ONE pass over the gzip stream -- the full decompressed
    # file never touches disk.
    loci_file = None
    try:
        with gzip.open(source_gz, "rt", newline="") as src, \
                open(normal_dest, "w", newline="\n") as normal_out, \
                open(tumor_dest, "w", newline="\n") as tumor_out:
            reader = csv.reader(src, delimiter="\t")
            normal_writer = csv.writer(normal_out, delimiter="\t", lineterminator="\n")
            tumor_writer = csv.writer(tumor_out, delimiter="\t", lineterminator="\n")

            header = next(reader)
            index = {name: i for i, name in enumerate(header)}
            normal_plan = column_plan(index, config.LOCUS_COLUMNS + histogram_columns("NORMAL_"))
            tumor_plan = column_plan(index, config.LOCUS_COLUMNS + histogram_columns("TUMOR_"))
            normal_writer.writerow([name for name, _ in normal_plan])
            tumor_writer.writerow([name for name, _ in tumor_plan])

            add_germline = GERMLINE_COLUMN not in index
            if add_germline:
                print(f"  {os.path.basename(source_gz)} has no {GERMLINE_COLUMN}; "
                      f"reading it line-for-line from {os.path.basename(config.GERMLINE_LOCI_FILE)}")
                loci_file = open(config.GERMLINE_LOCI_FILE)
                key_indices = [index[c] for c in config.LOCUS_KEY_COLUMNS]

            germline = "False"
            germline_count = rows = 0
            for row in reader:
                rows += 1
                if add_germline:
                    loci_line = loci_file.readline()
                    if not loci_line:
                        raise RuntimeError(
                            f"{config.GERMLINE_LOCI_FILE} ran out of loci at row {rows} of "
                            f"{source_gz}; the two files do not line up.")
                    loci_locus, annotation = parse_loci_line(loci_line)
                    if annotation not in VALID_LOCI_GERMLINE_VALUES:
                        raise RuntimeError(
                            f"unexpected germline annotation {annotation!r} at line {rows} of "
                            f"{config.GERMLINE_LOCI_FILE}; expected one of {VALID_LOCI_GERMLINE_VALUES}.")
                    mut_locus = tuple(row[i] for i in key_indices)
                    if mut_locus != loci_locus:
                        raise RuntimeError(
                            f"{source_gz} and the loci file disagree at row {rows}: "
                            f"mutation file has {mut_locus}, loci file has {loci_locus}. "
                            f"The line-for-line germline annotation is only valid when a "
                            f"*.full.mut.tsv has one row per locus in loci-file order.")
                    germline = "True" if annotation == HETEROGENEOUS else "False"
                    germline_count += annotation == HETEROGENEOUS
                normal_writer.writerow([germline if i is None else row[i] for _, i in normal_plan])
                tumor_writer.writerow([germline if i is None else row[i] for _, i in tumor_plan])
    finally:
        if loci_file is not None:
            loci_file.close()

    if add_germline:
        print(f"  {germline_count}/{rows} loci ({100 * germline_count / max(rows, 1):.1f}%) "
              f"marked {GERMLINE_COLUMN}=True")


def run_sample(sample: str, force: bool = False, keep_full_gz: bool = True) -> str:
    config.ensure_free_disk_gb()

    full_gz = config.full_gz_local_path(sample)
    if not os.path.exists(full_gz):
        download_samples.download_full(sample)

    out_dir = config.rerun_sample_dir(sample)
    os.makedirs(out_dir, exist_ok=True)
    output_prefix = config.rerun_output_prefix(sample)
    output_file = output_prefix + ".full.mut.tsv"
    if os.path.exists(output_file) and not force:
        print(f"[{sample}] already re-run, skipping (use --force to redo): {output_file}")
        return output_file

    normal_hist = os.path.join(out_dir, f"{sample}.pseudo_normal.hist.tsv")
    tumor_hist = os.path.join(out_dir, f"{sample}.pseudo_tumor.hist.tsv")
    print(f"[{sample}] carving histograms from {full_gz}")
    carve_histograms(full_gz, normal_hist, tumor_hist)

    # the actual msmutect invocation is layer 2's job
    print(f"[{sample}] running msmutect --from_file")
    try:
        run_single_sample.run_from_histograms(normal_hist, tumor_hist, output_prefix,
                                              msmutect=config.MSMUTECT_PATH)
    finally:
        for path in (normal_hist, tumor_hist):
            if os.path.exists(path):
                os.remove(path)

    if not keep_full_gz:
        os.remove(full_gz)

    print(f"[{sample}] done -> {output_file}")
    return output_file


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", nargs="+", default=config.SAMPLES)
    parser.add_argument("--force", action="store_true", help="re-run even if output already exists")
    parser.add_argument("--delete-full-gz", action="store_true",
                         help="delete the downloaded *.full.mut.tsv.gz after processing each sample "
                              "to save disk space (off by default so re-runs don't need re-downloading)")
    args = parser.parse_args()

    for sample in args.samples:
        run_sample(sample, force=args.force, keep_full_gz=not args.delete_full_gz)


if __name__ == "__main__":
    main()
