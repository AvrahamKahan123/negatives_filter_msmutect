import pandas as pd
import numpy as np
import os
import re
from glob import glob
from functools import reduce

def output_edit(df: pd.DataFrame, case_name: str):
    # Extracting called mutations only
    df = df[df["CALL"] == "M"].copy()
    df.fillna(0, inplace=True)
    df["Motif_size"] = df["PATTERN"].astype(str).str.len()
    df["Case"] = case_name
    return df


def allele_num(df, new_col, allele1, allele2=None, allele3=None, allele4=None):
    if allele1 == "TUMOR_ALLELE_1":
        df[new_col] = (
            (df[allele1] > 0).astype(int) +
            (df[allele2] > 0).astype(int) +
            (df[allele3] > 0).astype(int) +
            (df[allele4] > 0).astype(int)
        )
    elif allele1 == "NORMAL_ALLELE_1":
        df[new_col] = (
            (df[allele1] > 0).astype(int) +
            (df[allele2] > 0).astype(int)
        )
    return df


def BTR_remove(df):
    df["REFERENCE_REPEATS"] = np.floor(df["REFERENCE_REPEATS"])

    # Convert relevant columns to Series for cleaner code
    R = np.floor(df["REFERENCE_REPEATS"])
    N1 = df["NORMAL_ALLELE_1"]
    N2 = df["NORMAL_ALLELES_2"]
    T1 = df["TUMOR_ALLELE_1"]
    T2 = df["TUMOR_ALLELES_2"]
    T3 = df["TUMOR_ALLELES_3"]
    T4 = df["TUMOR_ALLELES_4"]

    # Condition 1: 1 normal allele
    cond1 = (
        (N2 == 0) &
        (
            ((T2 == 0) & (T1 == R)) |
            ((T3 == 0) & (T2 != 0) & (N1 != R) &
             ((T2 == R) | (T1 == R)) &
             ((T2 == N1) | (T1 == N1)))
        )
    )

    # Condition 2: 2 normal alleles and 2 tumor alleles
    cond2 = (
        (N2 != 0) & (T3 == 0) &
        (
            ((T1 == R) & (T2.isin([N1, N2]))) |
            ((T2 == R) & (T1.isin([N1, N2])))
        )
    )

    # Condition 3: 2 normal alleles and 3 tumor alleles
    cond3 = (
        (N2 != 0) & (T4 == 0) & (T3 != 0) &
        (
            ((T1 == R) & T2.isin([N1, N2]) & T3.isin([N1, N2])) |
            ((T2 == R) & T1.isin([N1, N2]) & T3.isin([N1, N2])) |
            ((T3 == R) & T1.isin([N1, N2]) & T2.isin([N1, N2]))
        )
    )

    remove_mask = cond1 | cond2 | cond3
    print(cond1.sum())
    print(cond2.sum())
    print(cond3.sum())

    # Return dataframe with BTR mutations removed

    return df[~remove_mask].copy()


def output_filt(df):
    # Removing copy number variations (loss of heterozygosity)
    df = allele_num(df, "norm_allele_num", "NORMAL_ALLELE_1", "NORMAL_ALLELES_2")

    df = allele_num(
        df,
        "tum_allele_num",
        "TUMOR_ALLELE_1",
        "TUMOR_ALLELES_2",
        "TUMOR_ALLELES_3",
        "TUMOR_ALLELES_4"
    )

    # Keep rows where tumor allele count is at least as high as normal allele count
    df = df[df["norm_allele_num"] <= df["tum_allele_num"]]

    # Remove back-to-reference (BTR) mutations
    df = BTR_remove(df)

    return df


def total_mut_count(df, case_name):
    return pd.DataFrame([[case_name, len(df)]], columns=["Case", "Total_mut_count"])


def per_motif_mut_count(df, case_name):
    row = {"Case": case_name}
    for n in range(1, 16):
        row[f"Motif_size{n}_mut_count"] = len(df[df["Motif_size"] == n])
    return pd.DataFrame([row])


def file_write(df, case_name):
    output_dir = "data/muts_only"
    os.makedirs(output_dir, exist_ok=True)

    output_path = os.path.join(output_dir, f"{case_name}.called.filt.mut.tsv.gz")
    df.to_csv(output_path, sep="\t", index=False, compression="gzip", quoting=3)  # quoting=3 means no quotes


def add_dfs(dfs):
    first_df = dfs[0]
    for df in dfs[1:]:
        first_df.iloc[0, 1:] = df.iloc[0, 1:] + first_df.iloc[0, 1:]
    return first_df


def process_individual_file(fp, case_name):
    total_mut_count_list = []
    per_motif_mut_count_list = []
    chunk_iter = pd.read_csv(fp, sep="\t", compression="gzip", chunksize=100_000)
    # chunk = next(chunk_iter)
    for chunk in chunk_iter:
    # for i in range(2):
        output_table = chunk
        output_table = output_edit(output_table, case_name)
        output_table = output_filt(output_table)

        file_write(output_table, case_name)
        total_mut_count_list.append(total_mut_count(output_table, case_name))
        per_motif_mut_count_list.append(per_motif_mut_count(output_table, case_name))
    return add_dfs(total_mut_count_list), add_dfs(per_motif_mut_count_list)
    # return reduce(lambda x, y: x.add(y, fill_value=0), total_mut_count_list), reduce(lambda x, y: x.add(y, fill_value=0), per_motif_mut_count_list)
    # return pd.concat(total_mut_count_list, ignore_index=True), pd.concat(per_motif_mut_count_list, ignore_index=True)


def create_final_table(input_dir):
    files = glob(os.path.join(input_dir, "**/*.tsv.gz"), recursive=True)

    case_names = [os.path.basename(f)[:-16] for f in files]

    total_mut_counts_list = []
    total_per_motif_mut_counts = []
    total_mut_count_table = pd.DataFrame(columns=["Case", "Total_mut_count"])
    motif_cols = [f"Motif_size{n}_mut_count" for n in range(1, 16)]
    per_motif_size_mut_count_table = pd.DataFrame(columns=["Case"] + motif_cols)

    # Step 4: Process each file
    for file, case_name in zip(files, case_names):
        mut_count, mut_count_per_motif = process_individual_file(file, case_name)
        total_mut_counts_list.append(mut_count)
        # total_mut_count_table = pd.concat(
        #     [total_mut_count_table, mut_count],
        #     ignore_index=True
        # )
        total_per_motif_mut_counts.append(mut_count_per_motif)
        # per_motif_size_mut_count_table = pd.concat(
        #     [per_motif_size_mut_count_table, mut_count_per_motif],
        #     ignore_index=True
        # )
    total_mut_count_table = pd.concat(total_mut_counts_list, ignore_index=True)
    per_motif_size_mut_count_table = pd.concat(total_per_motif_mut_counts, ignore_index=True)
    os.makedirs("Output_tables", exist_ok=True)

    total_mut_count_table.to_csv("Output_tables/total_mut_count_table.txt", sep="\t", index=False)
    per_motif_size_mut_count_table.to_csv("Output_tables/per_motif_size_mut_count_table.txt", sep="\t", index=False)


if __name__ == '__main__':
    ## CHANGE HERE ##
    create_final_table("Input_files_test")
    # pass
