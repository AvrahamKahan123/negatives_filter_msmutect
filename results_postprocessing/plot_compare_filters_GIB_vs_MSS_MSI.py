import os

import pandas as pd
import matplotlib
matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt
import seaborn as sns

# NEGATIVE_CONTROL rows are the GIB samples; MSS/MSI come from TCGA.
FILTER_COLUMNS = {"AVRAHAM_FILTER": "Avraham", "YOSSI_FILTER": "Yossi"}
GROUP_ORDER = ["GIB", "MSS", "MSI"]
FILTER_ORDER = ["Avraham", "Yossi"]
OUT_FP = "/home/avraham/MaruvkaLab/msmutect_postprocessing/results/compare_filters_GIB_vs_MSS_MSI.violin.png"


def load_long(mss_and_msi_source: str, gib_source: str) -> pd.DataFrame:
    df = pd.concat([pd.read_csv(gib_source), pd.read_csv(mss_and_msi_source)], ignore_index=True)
    df["GROUP"] = df["CLASSIFICATION"].replace({"NEGATIVE_CONTROL": "GIB"})
    # one row per (sample, filter): MUTATIONS = number of calls surviving that filter
    long = df.melt(id_vars=["CASE", "GROUP"], value_vars=list(FILTER_COLUMNS),
                   var_name="FILTER", value_name="MUTATIONS")
    long["FILTER"] = long["FILTER"].map(FILTER_COLUMNS)
    return long


def main():
    mss_and_msi_source = "/home/avraham/MaruvkaLab/msmutect_postprocessing/results/full_tcga_COMPARE_FILTERS_9_2.csv"
    gib_source = "/home/avraham/MaruvkaLab/msmutect_postprocessing/results/gib_COMPARE_FILTERS_9_2.csv"

    long = load_long(mss_and_msi_source, gib_source)

    # sample count per group (identical for both filters) -> shown under each group label
    n_per_group = long[long["FILTER"] == FILTER_ORDER[0]].groupby("GROUP")["CASE"].count()

    fig, ax = plt.subplots(figsize=(9, 6))
    sns.violinplot(
        data=long, x="GROUP", y="MUTATIONS", hue="FILTER",
        order=GROUP_ORDER, hue_order=FILTER_ORDER,
        log_scale=True,          # KDE computed in log space, y-axis log
        density_norm="width",    # equal max width per violin (group sizes differ a lot: 50 vs 7100 vs 286)
        cut=0, inner="quartile", ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel("mutations passing filter (log scale)")
    ax.set_title("Filter comparison: Avraham vs Yossi across GIB / MSS / MSI")
    ax.set_xticks(range(len(GROUP_ORDER)))
    ax.set_xticklabels([f"{g}\n(n={int(n_per_group.get(g, 0))})" for g in GROUP_ORDER])
    ax.legend(title="filter")
    fig.tight_layout()

    fig.savefig(OUT_FP, dpi=150)
    plt.close(fig)
    print("wrote", OUT_FP)


if __name__ == "__main__":
    main()
