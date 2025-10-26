import pandas as pd
import numpy as np


def is_deletion(row):
    norm1 = row['NORMAL_ALLELE_1']
    norm2 = row['NORMAL_ALLELES_2']
    tumor_alleles = [row.get(col) for col in ['TUMOR_ALLELE_1', 'TUMOR_ALLELE_2', 'TUMOR_ALLELE_3'] if pd.notna(row.get(col))]

    # Case a: norm2 is N/A or missing
    if pd.isna(norm2):
        for allele in tumor_alleles:
            if allele != norm1:
                return len(allele) < len(norm1)
        return False

    # Case b: both norm1 and norm2 have values
    for allele in tumor_alleles:
        if allele != norm1 and allele != norm2:
            return len(allele) < len(norm1) and len(allele) < len(norm2)
    return False


def get_matching_rows(file_path):
    df = pd.read_csv(file_path, sep='\t', dtype=str)

    # Standardize N/A values to np.nan
    df.replace("N/A", np.nan, inplace=True)

    # Filter for call == 'M' and pattern == 'A'
    filtered = df[(df['CALL'] == 'M') & (df['PATTERN'] == 'A')].copy()

    # Keep only rows where is_deletion is True
    matching_rows = filtered[filtered.apply(is_deletion, axis=1)]

    return matching_rows


# Example usage:
result_df = get_matching_rows("/data/muts_only/tst_all.tsv")
croc=1
# print(result_df)
# result_df.to_csv("matching_rows.tsv", sep='\t', index=False)
