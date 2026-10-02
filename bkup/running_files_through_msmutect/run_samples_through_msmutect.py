import os, csv, glob, subprocess
import time
from typing import List

# Re-runs MSMuTect on GIAB samples "from file": each input is a .full.mut.tsv,
# which already contains the per-locus normal and tumor histograms. We slice a
# pseudo-normal and pseudo-tumor histogram back out of it (same cuts as the
# MSMuTect from_file E2E test) and feed them to msmutect --from_file, which
# recomputes the calls. For same-sample GIAB pairs this measures the mutation
# (false-positive) rate.
#
#     cut -f1-19       IN > <name>.pseudo_normal.hist.tsv
#     cut -f1-6,29-41  IN > <name>.pseudo_tumor.hist.tsv
#     msmutect.sh -N <normal> -T <tumor> -m -A -O <name> --from_file

MSMUTECT_SH = "/home/avraham/MaruvkaLab/MSMuTect/msmutect.sh"

# 1-indexed columns of a .full.mut.tsv that make up each pseudo-histogram.
# (The filtered GIAB files carry 4 extra trailing columns, but the normal/tumor
# histogram blocks stay at these positions, so the cuts are unchanged.)
NORMAL_HIST_COLUMNS: List[int] = list(range(1, 20))                       # 1-19
TUMOR_HIST_COLUMNS: List[int] = list(range(1, 7)) + list(range(29, 42))  # 1-6, 29-41


def cut_columns(source: str, columns: List[int], destination: str):
    # cross-platform `cut -f<columns>`
    with open(source, newline="") as src, open(destination, "w", newline="\n") as dst:
        reader = csv.reader(src, delimiter="\t")
        writer = csv.writer(dst, delimiter="\t", lineterminator="\n")
        for row in reader:
            writer.writerow([row[column - 1] for column in columns])


def sample_name(fp: str) -> str:
    name = os.path.basename(fp)
    for suffix in (".filtered.full.mut.tsv", ".full.mut.tsv"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return os.path.splitext(name)[0]


def run_sample_through_msmutect(fp: str,
                                outdir="/home/avraham/MaruvkaLab/msmutect_postprocessing/results/filtered_through_msmutect_directly") -> str:
    os.makedirs(outdir, exist_ok=True)
    name = sample_name(fp)

    # 1. carve the pseudo-normal / pseudo-tumor histograms out of the one input file
    normal_hist = os.path.join(outdir, f"{name}.pseudo_normal.hist.tsv")
    tumor_hist = os.path.join(outdir, f"{name}.pseudo_tumor.hist.tsv")
    cut_columns(fp, NORMAL_HIST_COLUMNS, normal_hist)
    cut_columns(fp, TUMOR_HIST_COLUMNS, tumor_hist)

    # 2. recompute calls from those histograms with msmutect --from_file
    output_prefix = os.path.join(outdir, name)
    command = [
        MSMUTECT_SH,
        "-N", normal_hist,
        "-T", tumor_hist,
        "-m", "-A",
        "-O", output_prefix,
        "--from_file",
        "-f",  # overwrite on re-run
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


def run_dir_through_msmutect():
    dir_path = "/data/gib_files_filtered"
    pattern = "*.filtered.full.mut.tsv"
    for fp in glob.glob(os.path.join(dir_path, pattern)):
        st = time.time()
        run_sample_through_msmutect(fp)
        e=time.time()
        print(f"Ran: {os.path.basename(fp)} in {e-st}")


if __name__ == "__main__":
    # all
    run_dir_through_msmutect()

    # single test case
    # test_file = "/home/avraham/MaruvkaLab/msmutect_postprocessing/data/gib_files_filtered/hg001_normal0_tumor0.filtered.full.mut.tsv"
    # output_file = run_sample_through_msmutect(test_file)
    # print("wrote:", output_file)
