import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from enums import COLUMN, MSI_CLASSIFICATION


def main():
    gib_df = pd.read_csv("../results/gib_purity.csv")
    # gib_df = pd.read_csv("results/full_tcga_purity.csv")
    # gib_df = gib_df[gib_df[COLUMN.CLASSIFICATION]==MSI_CLASSIFICATION.MSS]
    min_size = 85
    dataset_name = "GIB"
    stat_name = "Purity"
    yticks = list(np.arange(0, 1.1, 0.1))
    # yticks = list(np.arange(0, 0.11, 0.01))


    per_patient_muts = gib_df.iloc[:, 3:]
    num_columns = len(per_patient_muts.columns)
    # Compute sums
    sums = per_patient_muts.sum()
    percentages = sums/sums.sum()
    plt.figure(figsize=(12, 6))
    percentages.plot(kind="bar", edgecolor="black", alpha=0.7)
    plt.title(f"{dataset_name}: {stat_name}")
    plt.xlabel(f"{stat_name}")
    plt.ylabel("Mutation Percentage")
    positions = list(range(0, len(per_patient_muts.columns)))
    plt.yticks(yticks)
    # labels = [str(i) for i in range(min_size, min_size+num_columns-1)] + ["100<"]
    labels = [str(i) for i in range(min_size, min_size+num_columns)]

    plt.xticks(positions, labels, rotation=90)
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    main()