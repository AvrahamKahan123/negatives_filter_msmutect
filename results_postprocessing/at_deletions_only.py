import os

import pandas
import pandas as pd
import numpy as np
from pathlib import Path


def ppdf(df: pandas.DataFrame):
    # pretty print data frame
    df.to_csv("tmp.tsv", sep="\t", index=False)
    print(Path("tmp.tsv").read_text())


def hom_to_het_extract(df, old=False):
    if old:
        NORMAL_ALLELE_2 = "NORMAL_ALLELES_2"
        TUMOR_ALLELE_3 = "TUMOR_ALLELES_3"
        TUMOR_ALLELE_2 = "TUMOR_ALLELES_2"
    else:
        NORMAL_ALLELE_2 = "NORMAL_ALLELE_2"
        TUMOR_ALLELE_3 = "TUMOR_ALLELE_3"
        TUMOR_ALLELE_2 = "TUMOR_ALLELE_2"


    # Filter rows matching the condition
    condition = (
        (df[NORMAL_ALLELE_2] == 0) &
        (df[TUMOR_ALLELE_3] == 0) &
        (df[TUMOR_ALLELE_2] != 0) &
        (
            (df["TUMOR_ALLELE_1"] == df["NORMAL_ALLELE_1"]) |
            (df[TUMOR_ALLELE_2] == df["NORMAL_ALLELE_1"])
        )
    )
    df = df[condition].copy()

    # Calculate Indel_length
    df["Indel_length"] = np.where(
        df["NORMAL_ALLELE_1"] == df["TUMOR_ALLELE_1"],
        df[TUMOR_ALLELE_2] - df["NORMAL_ALLELE_1"],
        df["TUMOR_ALLELE_1"] - df["NORMAL_ALLELE_1"]
    )

    return df[df["Indel_length"] < 0]



def process_directory(filtered_mut_files_dir: str, output_table_filename: str):
    """
    :param filtered_mut_files_dir: the directory holding all of the FILTERED .filtered.full.mut.tsv
    :param output_table_filename: the name of the outputted table
    :return:
    """
    results = []
    for filename in os.listdir(filtered_mut_files_dir):
        if filename.endswith(".mut.tsv"):
            file_path = os.path.join(filtered_mut_files_dir, filename)
            case_name = filename[:-len(".filtered.full.mut.tsv")]
            df = pd.read_csv(file_path, sep="\t")
            relevant_df = df[(df["PATTERN"].isin(["A", "T"])) & (df["CALL"] == "M")]
            filtered_df = hom_to_het_extract(relevant_df, old=True)
            results.append((case_name, len(filtered_df)))

    with open(output_table_filename, "w") as f:
        f.write("CASE_NAME\tA_T_DELETIONS\n")
        for case_name, count in results:
            f.write(f"{case_name}\t{count}\n")


if __name__ == '__main__':
    ## CHANGE HERE ##
    process_directory("../positives_simulation/output_files", "../positives_simulation/complex_simulation.txt")
