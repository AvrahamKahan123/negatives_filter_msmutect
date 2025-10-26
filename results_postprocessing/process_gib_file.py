import pandas as pd
import random

# # Example: load your dataframe
df = pd.read_csv("/home/avraham/MaruvkaLab/tmp/msmutect_normal0_tumor0.full.mut.tsv", sep="\t", nrows=100_000)
df = df.fillna(0)
df = df[df["REFERENCE_REPEATS"] % 1 == 0]
# Find rows where either NORMAL_ALLELE_1 or NORMAL_ALLELE_2 does not match REFERENCE_REPEATS
mask = (df["NORMAL_ALLELE_1"] != df["REFERENCE_REPEATS"]) & (df["NORMAL_ALLELE_2"] != df["REFERENCE_REPEATS"]) & (df["NORMAL_ALLELE_1"] !=0)

# Filter the dataframe
mismatched = df[mask]

print(len(mismatched))


