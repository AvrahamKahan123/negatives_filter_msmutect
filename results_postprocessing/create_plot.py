import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import re

# Read tables
total_mut_count_table = pd.read_csv("Output_tables/total_mut_count_table.txt", sep="\t")
per_motif_size_mut_count_table = pd.read_csv("Output_tables/per_motif_size_mut_count_table.txt", sep="\t")

# Melt the wide-format motif table to long-format
long_df = per_motif_size_mut_count_table.melt(
    id_vars="Case",
    var_name="Motif_size",
    value_name="Indel_count"
)

# Clean the Motif_size column: remove "Motif_size" and "_mut_count"
long_df["Motif_size"] = long_df["Motif_size"].str.replace("Motif_size", "", regex=False)
long_df["Motif_size"] = long_df["Motif_size"].str.replace("_mut_count", "", regex=False)
long_df["Motif_size"] = long_df["Motif_size"].astype(int)

# Filter to only motif sizes 1 through 8
motif_sizes_to_plot = list(range(1, 9))
plot_df = long_df[long_df["Motif_size"].isin(motif_sizes_to_plot)].copy()

# Add log10(Indel_count + 1) to avoid log(0)
plot_df["log10_Indel_count"] = np.log10(plot_df["Indel_count"] + 1)

# Plotting
plt.figure(figsize=(10, 6))
sns.violinplot(data=plot_df, x="Motif_size", y="log10_Indel_count", inner=None)
sns.boxplot(data=plot_df, x="Motif_size", y="log10_Indel_count", whis=1.5)
sns.stripplot(data=plot_df, x="Motif_size", y="log10_Indel_count", color='black', alpha=0.5)

# Styling to mimic theme_minimal + customizations
plt.ylabel("log10(Indel count)")
plt.xlabel("Motif size")
sns.despine()
plt.grid(False)
plt.xticks(color="black")
plt.yticks(color="black")
plt.tight_layout()
plt.show()