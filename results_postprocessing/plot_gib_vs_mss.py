from typing import Tuple
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from results_postprocessing.cancer_data import results_directory
from results_postprocessing.enums import COLUMN, MSI_CLASSIFICATION


def split_into_hg005_and_nonhg005(gib_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    gib_case_names = list(gib_df["CASE"])
    hg005_idxs = ["HG005" in case_name for case_name in gib_case_names]
    return gib_df[hg005_idxs], gib_df[[not(idx) for idx in hg005_idxs]]


def yossi_vs_avr_filter():
    gib_df: pd.DataFrame = pd.read_csv(os.path.join(results_directory(), "gib_COMPARE_FILTERS.csv"))
    _, gib_df = split_into_hg005_and_nonhg005(gib_df)
    df_melted = gib_df.melt(value_vars=['ORIGINAL', 'AVRAHAM_FILTER', 'YOSSI_FILTER'], var_name='Filter', value_name='Num Mutations')
    plt.yticks(list(range(0, 2000, 100)))
    # Plot
    sns.violinplot(x='Filter', y='Num Mutations', data=df_melted)
    plt.title("GIB response to different filters")
    plt.show()

    # all_gib_df = pd.read_csv("../results/gib_COMPARE_FILTERS.csv")
    # hg005_df, gib_df = split_into_hg005_and_nonhg005(all_gib_df)
    # gib_df[COLUMN.CLASSIFICATION] = "GIB"
    # # hg005_df[COLUMN.CLASSIFICATION] = "HG005"
    # combined_df = pd.concat([hg005_df, gib_df])
    # combined_df["LOG10_muts"] = np.log10(combined_df["13"])
    #
    #
    # # Make violin plot
    # sns.violinplot(x=COLUMN.CLASSIFICATION, y="LOG10_muts", data=combined_df)
    # plt.ylabel("Log10 number of mutations")
    # plt.yticks(list(range(1, 7)))
    #
    # plt.title(f"GIB vs MSS")
    #
    # plt.show()


def plot_mss_vs_gib(tcga_csv: str, gib_csv: str, column_name: str, title: str = "GIB vs MSS", save: bool = False):
    full_tcga_df = pd.read_csv(tcga_csv)
    mss_df = full_tcga_df[full_tcga_df[COLUMN.CLASSIFICATION]==MSI_CLASSIFICATION.MSS]
    all_gib_df = pd.read_csv(gib_csv)
    hg005_df, gib_df = split_into_hg005_and_nonhg005(all_gib_df)
    gib_df[COLUMN.CLASSIFICATION] = "GIB"
    hg005_df[COLUMN.CLASSIFICATION] = "HG005"
    combined_df = pd.concat([hg005_df, gib_df, mss_df])
    combined_df["LOG10_muts"] = np.log10(combined_df[column_name])


    # Make violin plot
    sns.violinplot(x=COLUMN.CLASSIFICATION, y="LOG10_muts", data=combined_df)
    plt.ylabel("Log10 number of mutations")
    plt.yticks(list(range(1, 7)))

    plt.title(title)
    if save:
        plt.savefig(f"{title}.png")
        plt.close()
    else:
        plt.show()


if __name__ == '__main__':
    for column in ["ORIGINAL"]+list(range(8, 16)):
        plot_mss_vs_gib("../results/full_tcga_NOISELESS_FILTER_NO_UTF.csv", "../results/gib_NOISELESS_FILTER_NO_UTF.csv",
                        f"{column}", f"{column}", save=True)
    # yossi_vs_avr_filter()