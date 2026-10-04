import csv, os, glob
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, ListedColormap
from matplotlib.patches import Patch

NODATA_COLOR = "#f0f0f0"  # no callable loci of this (pattern_length, num_repeats)
ZERO_COLOR = "#969696"    # callable, but zero mutations (a measured zero, not missing data)

MAX_PATTERN_LEN = 15
#MUTATED_FILES_DIR = "/home/avraham/MaruvkaLab/msmutect_postprocessing/results/filtered_through_msmutect_directly"
# MUTATED_FILES_DIR = "/home/avraham/MaruvkaLab/msmutect_postprocessing/results/gib_filtered_through_newest_msmutect"
MUTATED_FILES_DIR="/home/avraham/MaruvkaLab/msmutect_postprocessing/Negatives/count_GIAB_mutation_rates/results/with_fisher_fixed"

MUTATED_PATTERN = "*.full.mutated_only.mut.tsv"
COUNTED_DIR = "/home/avraham/MaruvkaLab/msmutect_postprocessing/results/callable_counts"
HEATMAP_DIR = "/home/avraham/MaruvkaLab/msmutect_postprocessing/results/mutation_density_heatmaps_10_2"


def to_int(value) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 0


def count_mutations_per_cell(mutated_msmutect_file: str) -> dict:
    # numerator: number of mutations per (pattern_length, num_repeats)
    counts = defaultdict(int)
    with open(mutated_msmutect_file, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row.get("CALL", "M").strip() != "M":
                continue
            pattern_len = len(row["PATTERN"].strip())
            num_repeats = to_int(row["REFERENCE_REPEATS"])
            counts[(pattern_len, num_repeats)] += 1
    return counts


def read_callable_per_cell(counted_callable_loci_file: str) -> dict:
    # CALLABLE (number of callable loci) per (pattern_length, num_repeats)
    callable_loci = {}
    with open(counted_callable_loci_file, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            key = (to_int(row["PATTERN_LEN"]), to_int(row["NUM_REPEATS"]))
            callable_loci[key] = to_int(row["CALLABLE"])
    return callable_loci


def sample_name(mutated_msmutect_file: str) -> str:
    return os.path.basename(mutated_msmutect_file).replace(".full.mutated_only.mut.tsv", "")


def heat_map_for_single_file(mutated_msmutect_file: str, counted_callable_loci_file: str) -> np.ndarray:
    """Mutation density heat map for one sample.

    density(pattern_length, num_repeats) = mutations / (pattern_length * num_repeats * num_callable)
    i.e. mutations per callable base of that (motif length, repeat count) type.
    Saves a PNG and returns the density matrix (rows: pattern_length 1..15, cols: num_repeats 1..max).
    """
    mutations = count_mutations_per_cell(mutated_msmutect_file)
    callable_loci = read_callable_per_cell(counted_callable_loci_file)

    max_repeats = max((nr for _, nr in callable_loci), default=0)
    density = np.full((MAX_PATTERN_LEN, max_repeats), np.nan)  # NaN = no callable bases -> density undefined
    placed = 0
    for (pattern_len, num_repeats), num_callable in callable_loci.items():
        if not (1 <= pattern_len <= MAX_PATTERN_LEN) or num_repeats < 1 or num_callable <= 0:
            continue
        total_bases = pattern_len * num_repeats * num_callable  # total bases of this type
        muts = mutations.get((pattern_len, num_repeats), 0)
        density[pattern_len - 1, num_repeats - 1] = muts / total_bases
        placed += muts

    total_muts = sum(mutations.values())
    if placed < total_muts:
        print(f"  note: {total_muts - placed}/{total_muts} mutations fell on (pattern_length, num_repeats) "
              f"cells with no callable count and were not placed")

    _render_heatmap(density, sample_name(mutated_msmutect_file), max_repeats)
    return density


def _render_heatmap(density: np.ndarray, name: str, max_repeats: int):
    os.makedirs(HEATMAP_DIR, exist_ok=True)
    fig, ax = plt.subplots(figsize=(max(6, max_repeats * 0.3), 6))

    # Layer 1 (background): distinguish "callable, 0 mutations" from "no callable loci".
    # Every callable cell (finite density, incl. exactly 0) is painted ZERO_COLOR; cells with
    # no callable loci (NaN) fall through to NODATA_COLOR.
    callable_bg = np.ma.masked_invalid(np.where(np.isfinite(density), 0.0, np.nan))
    zero_cmap = ListedColormap([ZERO_COLOR])
    zero_cmap.set_bad(NODATA_COLOR)
    ax.imshow(callable_bg, aspect="auto", origin="lower", cmap=zero_cmap, vmin=0, vmax=1)

    # Layer 2 (overlay): positive densities on a log scale; non-positive cells are transparent
    # so the background layer shows through.
    positive = np.ma.masked_where(~(density > 0), density)
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad(alpha=0.0)
    norm = None
    if positive.count():
        vmin, vmax = positive.min(), positive.max()
        if vmin == vmax:
            vmin = vmax / 10.0
        norm = LogNorm(vmin=vmin, vmax=vmax)
    im = ax.imshow(positive, aspect="auto", origin="lower", cmap=cmap, norm=norm)

    ax.set_xlabel("number of repeats")
    ax.set_ylabel("pattern length")
    ax.set_yticks(range(MAX_PATTERN_LEN))
    ax.set_yticklabels(range(1, MAX_PATTERN_LEN + 1))
    step = max(1, max_repeats // 20)
    ax.set_xticks(range(0, max_repeats, step))
    ax.set_xticklabels(range(1, max_repeats + 1, step))
    ax.set_title(f"{name}\nmutation density (mutations / callable base)")
    if positive.count():
        fig.colorbar(im, ax=ax, label="mutations per callable base")
    ax.legend(handles=[Patch(facecolor=ZERO_COLOR, edgecolor="gray", label="callable, 0 mutations"),
                       Patch(facecolor=NODATA_COLOR, edgecolor="gray", label="no callable loci")],
              loc="upper right", framealpha=0.9, fontsize=8)
    fig.tight_layout()

    out_fp = os.path.join(HEATMAP_DIR, f"{name}.heatmap.png")
    fig.savefig(out_fp, dpi=150)
    plt.close(fig)
    print(f"  wrote {out_fp}")


def count_file_for(mutated_fp: str) -> str:
    # mutated: hg001_normal0_tumor0.full.mutated_only.mut.tsv
    # count:   HG001_msmutect_normal0_tumor0.count.txt
    key = sample_name(mutated_fp).lower()
    key_no_suffix = key[:key.find(".")]
    for count_fp in glob.glob(os.path.join(COUNTED_DIR, "*.count.txt")):
        filename = os.path.basename(count_fp)
        count_sample_name = filename[:filename.find(".")].lower()
        if count_sample_name == key_no_suffix:
            return count_fp
        # count_key = os.path.basename(count_fp).replace(".count.txt", "").replace("msmutect_", "").lower()
        # if count_key == key:
        #     return count_fp
    return ""


def run_on_all():
    for mutated_fp in sorted(glob.glob(os.path.join(MUTATED_FILES_DIR, MUTATED_PATTERN))):
        count_fp = count_file_for(mutated_fp)
        if not count_fp:
            print(f"{sample_name(mutated_fp)}: no matching count file, skipping")
            continue
        print(sample_name(mutated_fp))
        heat_map_for_single_file(mutated_fp, count_fp)


def main():
    run_on_all()


if __name__ == "__main__":
    main()
