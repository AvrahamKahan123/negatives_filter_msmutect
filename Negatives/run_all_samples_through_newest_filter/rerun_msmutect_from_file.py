import argparse
import csv
import gzip
import os
import subprocess

import config
import download_samples

# Step 3 (see Plan.md): run the newest MSMuTect on each sample's *.full.mut.tsv.gz via
# `--from_file`, the same way run_all_GIB_files_through_new_msmutect/run_msmutect.py does it:
# each input already holds the per-locus normal and tumor histograms, so we carve a
# pseudo-normal and pseudo-tumor histogram out of it (BY COLUMN NAME) and feed them to
# `msmutect --from_file`, which recomputes CALL/FISHER_TEST_P_VALUE/etc. with the current code.
#
# Samples are processed one at a time, and histograms are carved by streaming straight out
# of the downloaded .gz (never writing an uncompressed copy of the full file to disk), to
# keep peak disk usage low -- see Plan.md's disk-space warning.


def histogram_columns(prefix: str):
    return ([f"{prefix}MOTIF_REPEATS_{i}" for i in range(1, 7)] +
            [f"{prefix}SUPPORTING_READS_{i}" for i in range(1, 7)])


def cut_columns_from_gz(source_gz: str, wanted_columns, destination: str):
    # select the wanted columns BY NAME from the source header; skip any not present.
    # reads directly from the gzip stream so the full decompressed file never touches disk.
    with gzip.open(source_gz, "rt", newline="") as src, open(destination, "w", newline="\n") as dst:
        reader = csv.reader(src, delimiter="\t")
        writer = csv.writer(dst, delimiter="\t", lineterminator="\n")
        header = next(reader)
        index = {name: i for i, name in enumerate(header)}
        keep = [index[name] for name in wanted_columns if name in index]
        writer.writerow([header[i] for i in keep])
        for row in reader:
            writer.writerow([row[i] for i in keep])


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
    cut_columns_from_gz(full_gz, config.LOCUS_COLUMNS + histogram_columns("NORMAL_"), normal_hist)
    cut_columns_from_gz(full_gz, config.LOCUS_COLUMNS + histogram_columns("TUMOR_"), tumor_hist)

    command = [
        config.MSMUTECT_PATH,
        "-N", normal_hist,
        "-T", tumor_hist,
        "-m", "-A",
        "-O", output_prefix,
        "--from_file",
        "-f",
    ]
    print(f"[{sample}] running: {' '.join(command)}")
    result = subprocess.run(command, capture_output=True, text=True)
    os.remove(normal_hist)
    os.remove(tumor_hist)
    if result.returncode != 0:
        raise RuntimeError(f"msmutect from_file run failed ({result.returncode}) for {sample}\n"
                            f"command: {' '.join(command)}\nstderr:\n{result.stderr}")
    if not os.path.exists(output_file):
        raise RuntimeError(f"expected output not created: {output_file}")

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
