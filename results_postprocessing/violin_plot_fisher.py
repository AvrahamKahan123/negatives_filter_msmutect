import math

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np


def load_csv_without_zeros(fp: str):
    ret = pd.read_csv(fp, delimiter="\t")
    # ret["C"] = ret["C"].replace(0, 1)
    ret.iloc[:, 1:] = ret.iloc[:, 1:].replace(0, 1)
    return ret


def cles_bruteforce(x, y):
    """
    Compute Common Language Effect Size (CLES) using pairwise comparisons.
    """
    x = np.array(x)
    y = np.array(y)

    # Compare all pairs
    greater = np.sum(x[:, None] > y[None, :])
    equal = np.sum(x[:, None] == y[None, :])  # handle ties fairly
    total = x.size * y.size

    return (greater + 0.5 * equal) / total


# Load the CSVs
# gib = load_csv_without_zeros("gib_fisher.tsv")
# mss = load_csv_without_zeros("mss_fisher.tsv")
# msi = load_csv_without_zeros("msi_fisher.tsv")

gib = load_csv_without_zeros("results/gib_ks.tsv")
mss = load_csv_without_zeros("results/mss_ks.tsv")
msi = load_csv_without_zeros("results/msi_ks.tsv")

cols = ['ORIGINAL', 'LOG10 P_VAL=-0.25', 'LOG10 P_VAL=-0.5',
       'LOG10 P_VAL=-0.75', 'LOG10 P_VAL=-1.0', 'LOG10 P_VAL=-1.25',
       'LOG10 P_VAL=-1.5', 'LOG10 P_VAL=-1.75', 'LOG10 P_VAL=-2.0',
       'LOG10 P_VAL=-2.25', 'LOG10 P_VAL=-2.5', 'LOG10 P_VAL=-2.75',
       'LOG10 P_VAL=-3.0', 'LOG10 P_VAL=-3.25', 'LOG10 P_VAL=-3.5',
       'LOG10 P_VAL=-3.75', 'LOG10 P_VAL=-4.0', 'LOG10 P_VAL=-4.25',
       'LOG10 P_VAL=-4.5', 'LOG10 P_VAL=-4.75', 'LOG10 P_VAL=-5.0',
       'LOG10 P_VAL=-5.25', 'LOG10 P_VAL=-5.5', 'LOG10 P_VAL=-5.75',
       'LOG10 P_VAL=-6.0', 'LOG10 P_VAL=-6.25', 'LOG10 P_VAL=-6.5',
       'LOG10 P_VAL=-6.75', 'LOG10 P_VAL=-7.0']
for col in cols:

# Get column 2 (second column, index 1)
    gib_col = np.log10(gib[col])
    mss_col = np.log10(mss[col])
    msi_col = np.log10(msi[col])

    # Combine into one DataFrame with a label
    combined = pd.DataFrame({
        "LOG10 Mutations": pd.concat([gib_col, mss_col, msi_col], ignore_index=True),
        "Dataset": (["GIB"] * len(gib_col)
                  + ["MSS"] * len(mss_col)
                  + ["MSI"] * len(msi_col))
    })

    # Plot violin plot
    CLES = cles_bruteforce(msi_col, mss_col)
    plt.figure(figsize=(8, 6))
    sns.violinplot(x="Dataset", y="LOG10 Mutations", data=combined)
    plt.title("LOG10 Num Mutations")
    plt.yticks(list(np.arange(0, 7.5, 0.5)))
    plt.text(0.01, 0.99, f"CLES={-math.log10(1-CLES)}",
             transform=plt.gca().transAxes,
             ha="left", va="top",
             fontsize=12, color="black")

    # plt.savefig(f"violin_plot_ks_{col}.png")
    plt.show()
    plt.close()
