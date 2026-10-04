
# The goal is to plot the TCGA vs Genome in a bottle results:
# Genome in a bottle results: /home/avraham/MaruvkaLab/msmutect_postprocessing/results_postprocessing/data/M_counts_GIAB.tsv
# TCGA results: /home/avraham/MaruvkaLab/msmutect_postprocessing/results_postprocessing/data/M_counts_all_10_3_26.tsv
# First step: add a cancer type, and MSS classification to the TCGA file: use ./results/full_tcga_NOISELESS_FILTER_NO_UTF.csv and merge. Make a new file with this result
# Second step: plot a side by side triple violinplot: plot 1: GIAB plot, plot 2: MSS plot, plot 3: MSI plot
# Third step: For each cancer type plot a double violinplot of the MSS vs MSI for that cancer type

import os

import pandas as pd
import matplotlib
matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt
import seaborn as sns

from results_postprocessing.enums import COLUMN, MSI_CLASSIFICATION

DATA_DIR = "/home/avraham/MaruvkaLab/msmutect_postprocessing/results_postprocessing/data"
OUT_DIR = "/home/avraham/MaruvkaLab/msmutect_postprocessing/results_postprocessing/results"

GIAB_COUNTS_FP = os.path.join(DATA_DIR, "M_counts_GIAB.tsv")
TCGA_COUNTS_FP = os.path.join(DATA_DIR, "M_counts_all_10_3_26.tsv")
TCGA_METADATA_FP = "/home/avraham/MaruvkaLab/msmutect_postprocessing/results/full_tcga_NOISELESS_FILTER_NO_UTF.csv"
ANNOTATED_TCGA_FP = os.path.join(OUT_DIR, "M_counts_all_10_3_26.annotated.tsv")

# HG005 is excluded: it is the son of the HG006/HG007 trio and behaves differently from the
# other GIAB samples, so earlier comparisons (plot_gib_vs_mss.py) keep it out of the GIB group.
EXCLUDED_GIAB_SAMPLES = "HG005"

GIAB = "GIAB"
MUTATIONS = "M_count"
GROUP = "GROUP"
SAMPLE = "sample"

GROUP_ORDER = [GIAB, MSI_CLASSIFICATION.MSS, MSI_CLASSIFICATION.MSI]
CLASSIFICATION_ORDER = [MSI_CLASSIFICATION.MSS, MSI_CLASSIFICATION.MSI]


def annotate_tcga_counts() -> pd.DataFrame:
    """First step: attach CANCER_TYPE and MSS/MSI CLASSIFICATION to the per sample M counts."""
    counts_df = pd.read_csv(TCGA_COUNTS_FP, sep="\t")
    metadata_df = pd.read_csv(TCGA_METADATA_FP, usecols=[COLUMN.CASE, COLUMN.CANCER_TYPE, COLUMN.CLASSIFICATION])
    annotated_df = counts_df.merge(metadata_df, left_on=SAMPLE, right_on=COLUMN.CASE, how="left")
    annotated_df = annotated_df.drop(columns=[COLUMN.CASE])

    unannotated = annotated_df[COLUMN.CLASSIFICATION].isna()
    if unannotated.any():
        print(f"WARNING: {int(unannotated.sum())} samples missing from {os.path.basename(TCGA_METADATA_FP)}, dropping them")
        annotated_df = annotated_df[~unannotated]

    os.makedirs(OUT_DIR, exist_ok=True)
    annotated_df.to_csv(ANNOTATED_TCGA_FP, sep="\t", index=False)
    print("wrote", ANNOTATED_TCGA_FP)
    return annotated_df


def load_giab_counts() -> pd.DataFrame:
    giab_df = pd.read_csv(GIAB_COUNTS_FP, sep="\t")
    return giab_df[~giab_df[SAMPLE].str.contains(EXCLUDED_GIAB_SAMPLES)]


def violinplot(df: pd.DataFrame, x: str, order: list, title: str, out_fp: str):
    """One violin per value of x, mutation counts on a log y axis."""
    counts_per_group = df.groupby(x, observed=True)[SAMPLE].count()

    fig, ax = plt.subplots(figsize=(8, 6))
    sns.violinplot(
        data=df, x=x, y=MUTATIONS, hue=x, order=order, hue_order=order, legend=False,
        log_scale=True,          # KDE computed in log space, y axis log
        density_norm="width",    # equal max width per violin (group sizes differ by orders of magnitude)
        cut=0, inner="quartile", ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel("number of mutations (log scale)")
    ax.set_title(title)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([f"{group}\n(n={int(counts_per_group.get(group, 0))})" for group in order])
    fig.tight_layout()

    fig.savefig(out_fp, dpi=150)
    plt.close(fig)
    print("wrote", out_fp)


def plot_giab_vs_mss_vs_msi(giab_df: pd.DataFrame, tcga_df: pd.DataFrame):
    """Second step: GIAB vs MSS vs MSI, side by side."""
    giab_df = giab_df.assign(**{GROUP: GIAB})
    tcga_df = tcga_df.assign(**{GROUP: tcga_df[COLUMN.CLASSIFICATION]})
    combined_df = pd.concat([giab_df, tcga_df], ignore_index=True)
    violinplot(combined_df, GROUP, GROUP_ORDER, "GIAB vs TCGA MSS vs TCGA MSI",
               os.path.join(OUT_DIR, "GIAB_vs_MSS_vs_MSI.violin.png"))


def plot_mss_vs_msi_per_cancer_type(tcga_df: pd.DataFrame):
    """Third step: one MSS vs MSI figure per cancer type that has at least one MSI sample."""
    for cancer_type, cancer_type_df in tcga_df.groupby(COLUMN.CANCER_TYPE):
        n_msi = int((cancer_type_df[COLUMN.CLASSIFICATION] == MSI_CLASSIFICATION.MSI).sum())
        if n_msi == 0:
            print(f"skipping {cancer_type}: no MSI samples")
            continue
        violinplot(cancer_type_df, COLUMN.CLASSIFICATION, CLASSIFICATION_ORDER, f"{cancer_type}: MSS vs MSI",
                   os.path.join(OUT_DIR, f"MSS_vs_MSI.{cancer_type}.violin.png"))


def main():
    tcga_df = annotate_tcga_counts()
    giab_df = load_giab_counts()
    plot_giab_vs_mss_vs_msi(giab_df, tcga_df)
    plot_mss_vs_msi_per_cancer_type(tcga_df)


if __name__ == "__main__":
    main()
