"""Sensitivity of MSMuTect to tumor purity, for the 5 simulated MSI samples.

For each sample the simulated tumor histograms were rebuilt at purities 0.1 to 1.0
(see README.md) and re-called. This plots, per purity,

    true positive rate = mutations called in the simulated tumor
                         ---------------------------------------
                         mutations called in the real sample

The denominator comes from the whole-TCGA M counts, so numerator and denominator
are both current-MSMuTect calls on the same sample and the curve saturates near 1.

Writes a vector PDF for submission, a 600 dpi PNG, and the plotted numbers as CSV.
"""
import csv
import os
import re
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
data_file = "/home/avraham/MaruvkaLab/msmutect_postprocessing/positives_simulation/complex_simulation/results/M_counts_all.tsv"
denominator_file = "/home/avraham/MaruvkaLab/msmutect_postprocessing/results_postprocessing/data/M_counts_all_10_3_26.tsv"

OUT_PDF = os.path.join(HERE, "tpr_vs_purity_all_5_cases.pdf")
OUT_PNG = os.path.join(HERE, "tpr_vs_purity_all_5_cases.png")
OUT_CSV = os.path.join(HERE, "tpr_vs_purity_all_5_cases.csv")

# Validated against the data-viz six checks on the light surface #fcfcfb:
# every hue clears 3:1 contrast, worst adjacent CVD dE 13.3 (pass).
# Worst all-pair dE is 11.2, the floor band, so each series also carries a
# distinct marker -- colour is never the only thing separating two curves.
SERIES_COLORS = ["#008300", "#2a78d6", "#eb6834", "#4a3aa7", "#e34948"]
SERIES_MARKERS = ["o", "s", "^", "D", "v"]

TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"

# purity is the trailing part of e.g. TCGA-A6-5661_1x_0.35purity
NAME_RE = re.compile(r"^(?P<sample>.+?)_(?P<replicates>\d+)x_(?P<purity>[0-9.]+)purity$")


def read_counts(path):
    """sample id -> count, from a 2 column tsv with a header."""
    counts = {}
    with open(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader, None)
        for row in reader:
            if len(row) < 2 or not row[0].strip():
                continue
            counts[row[0].strip()] = int(row[1])
    return counts


def read_numerators(path):
    """sample -> {purity: count}, parsed out of the simulated file names."""
    curves = defaultdict(dict)
    unparsed = []
    for name, count in read_counts(path).items():
        match = NAME_RE.match(name)
        if match is None:
            unparsed.append(name)
            continue
        curves[match["sample"]][float(match["purity"])] = count
    if unparsed:
        sys.exit(f"could not parse purity out of: {unparsed}")
    return curves


def main():
    curves = read_numerators(data_file)
    totals = read_counts(denominator_file)

    # Sort by sample name, not by height, so a sample keeps its colour no matter
    # what the curves do or which subset is plotted.
    samples = sorted(curves)
    missing = [s for s in samples if s not in totals]
    if missing:
        sys.exit(f"no denominator in {denominator_file} for: {missing}")

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 9,
        "axes.linewidth": 0.8,
        "pdf.fonttype": 42,   # embed as editable TrueType, not outlines
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })

    # 180 mm, the usual double column width.
    fig, ax = plt.subplots(figsize=(7.09, 4.33))

    for i, sample in enumerate(samples):
        purities = sorted(curves[sample])
        rates = [100.0 * curves[sample][p] / totals[sample] for p in purities]
        ax.plot(
            purities, rates,
            color=SERIES_COLORS[i % len(SERIES_COLORS)],
            marker=SERIES_MARKERS[i % len(SERIES_MARKERS)],
            markersize=4.5,
            markeredgewidth=0,
            linewidth=1.8,
            label=f"{sample}  (n = {totals[sample]:,})",
            clip_on=False,
            zorder=3 + i,
        )

    ax.set_xlabel("Tumor purity", color=TEXT_PRIMARY)
    ax.set_ylabel("Mutations recovered (%)", color=TEXT_PRIMARY)
    ax.set_xlim(0.1, 1.0)
    ax.set_ylim(0, 100)
    ax.set_xticks([round(0.1 * i, 1) for i in range(1, 11)])
    ax.set_yticks(range(0, 101, 20))

    ax.yaxis.grid(True, color="#e6e5e1", linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(TEXT_SECONDARY)
    ax.tick_params(colors=TEXT_SECONDARY, labelcolor=TEXT_PRIMARY, length=3, width=0.8)

    legend = ax.legend(
        frameon=False,
        loc="upper left",
        handlelength=2.4,
        labelspacing=0.5,
        borderpad=0,
    )
    for text in legend.get_texts():
        text.set_color(TEXT_PRIMARY)

    fig.tight_layout()
    fig.savefig(OUT_PDF, bbox_inches="tight")
    fig.savefig(OUT_PNG, dpi=600, bbox_inches="tight")
    plt.close(fig)

    all_purities = sorted({p for s in samples for p in curves[s]})
    with open(OUT_CSV, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["purity"] + samples)
        for p in all_purities:
            row = [p]
            for s in samples:
                count = curves[s].get(p)
                row.append("" if count is None else round(100.0 * count / totals[s], 4))
            writer.writerow(row)

    print(f"{len(samples)} samples, {len(all_purities)} purities")
    for s in samples:
        top = max(curves[s])
        print(f"  {s}  denominator {totals[s]:>8,}  max {100.0 * curves[s][top] / totals[s]:5.1f}% at purity {top}")
    print(f"wrote {OUT_PDF}\n      {OUT_PNG}\n      {OUT_CSV}")


if __name__ == "__main__":
    main()
