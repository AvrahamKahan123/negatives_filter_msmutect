import os
import csv
import glob
import subprocess
import multiprocessing as mp

PARALLEL_RUNS = 6  # number of msmutect runs to execute concurrently

# Runs the newest MSMuTect on every new-format .full.mut.tsv in the GIB directory,
# "from file": each input already holds the per-locus normal and tumor histograms,
# so we carve a pseudo-normal and pseudo-tumor histogram out of it (BY COLUMN NAME)
# and feed them to `msmutect --from_file`, which recomputes the calls.
#
# The name-based carving + --from_file invocation is reused from
# MSMuTect/tests/E2E_tests/test_from_file_pre_called.py.

MSMUTECT_PATH = os.path.expanduser("~/MaruvkaLab/MSMuTect/msmutect.sh")
input_dir = "/home/avraham/MaruvkaLab/msmutect_postprocessing/data/gib_files_filtered_new_format"
input_pattern = "*.filtered.full.mut.tsv"
output_dir = \
    "/home/avraham/MaruvkaLab/msmutect_postprocessing/results/gib_filtered_through_newest_msmutect"

# locus columns include KNOWN_GERMLINE_VARIANT (the new from_file run requires it in
# the normal histogram); absent columns are simply skipped.
LOCUS_COLUMNS = ["CHROMOSOME", "START", "END", "PATTERN",
                 "REFERENCE_SEQUENCE", "REFERENCE_REPEATS", "KNOWN_GERMLINE_VARIANT"]


def histogram_columns(prefix: str):
    return ([f"{prefix}MOTIF_REPEATS_{i}" for i in range(1, 7)] +
            [f"{prefix}SUPPORTING_READS_{i}" for i in range(1, 7)])


def cut_columns(source: str, wanted_columns, destination: str):
    # select the wanted columns BY NAME from the source header; skip any not present
    with open(source, newline="") as src, open(destination, "w", newline="\n") as dst:
        reader = csv.reader(src, delimiter="\t")
        writer = csv.writer(dst, delimiter="\t", lineterminator="\n")
        header = next(reader)
        index = {name: i for i, name in enumerate(header)}
        keep = [index[name] for name in wanted_columns if name in index]
        writer.writerow([header[i] for i in keep])
        for row in reader:
            writer.writerow([row[i] for i in keep])


def sample_name(fp: str) -> str:
    name = os.path.basename(fp)
    for suffix in (".filtered.full.mut.tsv", ".full.mut.tsv"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return os.path.splitext(name)[0]


def run_sample_through_msmutect(fp: str, outdir: str = output_dir) -> str:
    os.makedirs(outdir, exist_ok=True)
    name = sample_name(fp)

    normal_hist = os.path.join(outdir, f"{name}.pseudo_normal.hist.tsv")
    tumor_hist = os.path.join(outdir, f"{name}.pseudo_tumor.hist.tsv")
    cut_columns(fp, LOCUS_COLUMNS + histogram_columns("NORMAL_"), normal_hist)
    cut_columns(fp, LOCUS_COLUMNS + histogram_columns("TUMOR_"), tumor_hist)

    output_prefix = os.path.join(outdir, name)
    command = [
        MSMUTECT_PATH,
        "-N", normal_hist,
        "-T", tumor_hist,
        "-m", "-A",
        "-O", output_prefix,
        "--from_file",
        "-f",
    ]
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"msmutect from_file run failed ({result.returncode}) for {fp}\n"
                           f"command: {' '.join(command)}\nstderr:\n{result.stderr}")
    output_file = output_prefix + ".full.mut.tsv"
    if not os.path.exists(output_file):
        raise RuntimeError(f"expected output not created: {output_file}")
    os.remove(normal_hist)
    os.remove(tumor_hist)
    return output_file


def run_all():
    files = sorted(glob.glob(os.path.join(input_dir, input_pattern)))
    os.makedirs(output_dir, exist_ok=True)
    # each --from_file run is single-core, so run PARALLEL_RUNS of them at once
    with mp.Pool(processes=PARALLEL_RUNS) as pool:
        for i, output_file in enumerate(pool.imap_unordered(run_sample_through_msmutect, files), 1):
            print(f"[{i}/{len(files)}] {os.path.basename(output_file)}")
    print(f"done: {len(files)} files -> {output_dir}")


if __name__ == "__main__":
    run_all()
