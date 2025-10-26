import os
from typing import List
import matplotlib.pyplot as plt
import numpy as np

import pandas as pd

import warnings
warnings.filterwarnings("ignore")

def allele_name(old=False):
    if old:
        return "ALLELES"
    else:
        return "ALLELE"


def second_third_normal_allele_counts(fps, title, old=True):
    print(f"*************\n{title}")
    for fp in fps:
        print(os.path.basename(fp))
        df = pd.read_csv(fp, delimiter="\t")
        num_muts = len(df)
        print(f"{len(df[df['NORMAL_MOTIF_REPEATS_2']!=0])/num_muts}, {len(df[df['NORMAL_MOTIF_REPEATS_3']!=0])/num_muts}")
        # print(f"NUMBER OF SECOND ALLELES: {len(df[df['NORMAL_MOTIF_REPEATS_2']!=0])}")
        # print(f"NUMBER OF THIRD ALLELES: {len(df[df['NORMAL_MOTIF_REPEATS_3']!=0])}")


def calc_median_allele_fraction(N1: int, N2: int, N3: int, N4: int, T1: int, T1S: int, T2: int, T2S: int, T3: int, T3S: int,
                                T4: int, T4S: int):
    normal_motifs = [N1, N2, N3, N4]
    tumor_motifs = [T1, T2, T3, T4]
    tumor_motif_support = [T1S, T2S, T3S, T4S]
    unique = 0
    total = sum(tumor_motif_support)
    for i in range(len(normal_motifs)):
        if tumor_motifs[i] not in normal_motifs:
            unique+=tumor_motif_support[i]
    return unique/total


def mutated_allele_fraction(fps: List[str], title: str) -> List[pd.DataFrame]:
    ret = []
    print(f"*************\n{title}")
    for fp in fps:
        print(os.path.basename(fp))
        df = pd.read_csv(fp, delimiter="\t")
        df["UNIQUE_TUMOR"] = df.apply(lambda row: calc_median_allele_fraction(row["NORMAL_MOTIF_REPEATS_1"],
                                                                              row["NORMAL_MOTIF_REPEATS_2"],
                                                                              row["NORMAL_MOTIF_REPEATS_3"],
                                                                              row["NORMAL_MOTIF_REPEATS_4"],

                                                                              row["TUMOR_MOTIF_REPEATS_1"],
                                                                              row["TUMOR_SUPPORTING_READS_1"],

                                                                              row["TUMOR_MOTIF_REPEATS_2"],
                                                                              row["TUMOR_SUPPORTING_READS_2"],

                                                                              row["TUMOR_MOTIF_REPEATS_3"],
                                                                              row["TUMOR_SUPPORTING_READS_3"],

                                                                              row["TUMOR_MOTIF_REPEATS_4"],
                                                                              row["TUMOR_SUPPORTING_READS_4"],

                                                                              ),
                                                                                axis=1)
        print(df["UNIQUE_TUMOR"].median())
        ret.append(df)
    return ret


def get_filename_no_ext(fp: str):
    fp = os.path.basename(fp)
    period_loc = fp.find(".")
    return fp[:period_loc]


def graph_column(fps: List[str], column_name: str, title: str, x_label: str, bins):
    for fp in fps:
        print(os.path.basename(fp))
        df = pd.read_csv(fp, delimiter="\t")
        data = df[column_name].dropna()
        data = data[(data != 0) & data.notna()]
        plt.figure(figsize=(8, 6))
        plt.hist(np.log10(data), bins=bins, edgecolor='black')
        plt.title(f"{title}_{get_filename_no_ext(fp)}")
        plt.xlabel(x_label)
        plt.ylabel("Num Mutations")

        # Save the figure
        plt.savefig(f"graphs/{title}_{get_filename_no_ext(fp)}.png")
        plt.close()


def graph_tumor_metric(fps: List[str], column_name: str, title: str, x_label: str, bins):
    dfs = mutated_allele_fraction(fps, "")
    for fp, df in zip(fps, dfs):
        print(os.path.basename(fp))
        # df = pd.read_csv(fp, delimiter="\t")
        data = df[column_name].dropna()
        data = data[(data != 0) & data.notna()]
        plt.figure(figsize=(8, 6))
        plt.hist(data, bins=bins, edgecolor='black')
        plt.title(f"{title}_{get_filename_no_ext(fp)}")
        plt.xlabel(x_label)
        plt.ylabel("Num Mutations")

        # Save the figure
        plt.savefig(f"graphs/{title}_{get_filename_no_ext(fp)}.png")
        plt.close()


def main():
    mss_files = [os.path.join("/data/mss", fp) for fp in
                 ['TCGA-2G-AALO.called.filt.mut.tsv',  'TCGA-D8-A27H.called.filt.mut.tsv',  'TCGA-E8-A414.called.filt.mut.tsv',
                  'TCGA-GM-A2DO.called.filt.mut.tsv',  'TCGA-LL-A5YO.called.filt.mut.tsv']]
    msi_files = [os.path.join("/data/msi", fp) for fp in
                 ['TCGA-A6-5661.called.filt.mut.tsv',  'TCGA-AJ-A3BH.called.filt.mut.tsv',  'TCGA-AP-A05N.called.filt.mut.tsv',
                  'TCGA-FI-A2D4.called.filt.mut.tsv',  'TCGA-OR-A5LB.called.filt.mut.tsv']]
    gib_files = ["C:/Users/avrah/MaruvkaLab/post_processing_code_for_australians/data/hg001_filtered/msmutect_normal2_tumor1.filtered.full.mut.tsv",
                 "C:/Users/avrah/MaruvkaLab/post_processing_code_for_australians/data/hg001_filtered/msmutect_normal0_tumor0.filtered.full.mut.tsv"]
    gib_files_2 = [
        "C:/Users/avrah/MaruvkaLab/post_processing_code_for_australians/data/hg002_filtered/msmutect_normal1_tumor0.filtered.full.mut.tsv",
        "C:/Users/avrah/MaruvkaLab/post_processing_code_for_australians/data/hg002_filtered/msmutect_normal2_tumor1.filtered.full.mut.tsv"]
    gib_files_3 = [
        "C:/Users/avrah/MaruvkaLab/post_processing_code_for_australians/data/hg003_filtered/msmutect_normal3_tumor1.filtered.full.mut.tsv",
        "C:/Users/avrah/MaruvkaLab/post_processing_code_for_australians/data/hg003_filtered/msmutect_normal2_tumor0.filtered.full.mut.tsv"]

    # second_third_normal_allele_counts(mss_files, "MSS")
    # second_third_normal_allele_counts(msi_files, "MSI")

    metric =["TUMOR_FRACTION"]#, "TUMOR_FRACTION", "FISHER"]
    #####################
    # KS TEST PVALUE
    if "KS" in metric:
        bins = np.arange(-20, 0, 0.25)

        graph_column(msi_files, "KS_TEST_PVALUE", "KS_TEST_PVALUE_msi", "KS_TEST_PVALUE_LOG10", bins)
        graph_column(mss_files, "KS_TEST_PVALUE", "KS_TEST_PVALUE_mss", "KS_TEST_PVALUE_LOG10", bins)
        graph_column(gib_files_2, "KS_TEST_PVALUE", "KS_TEST_PVALUE_HG002", "KS_TEST_PVALUE_LOG10", bins)
        graph_column(gib_files_3, "KS_TEST_PVALUE", "KS_TEST_PVALUE_HG003", "KS_TEST_PVALUE_LOG10", bins)

    if "FISHER" in metric:
        #FISHER TEST PVALUE
        bins = np.arange(-10, 0, 0.2)

        graph_column(msi_files, "FISHER_TEST_P_VALUE", "FISHER_PVALUE_msi", "FISHER_PVALUE_LOG10", bins)
        graph_column(mss_files, "FISHER_TEST_P_VALUE", "FISHER_PVALUE_mss", "FISHER_PVALUE_LOG10", bins)
        graph_column(gib_files_2, "FISHER_TEST_P_VALUE", "FISHER_PVALUE_HG002", "FISHER_PVALUE_LOG10", bins)
        graph_column(gib_files_3, "FISHER_TEST_P_VALUE", "FISHER_PVALUE_HG003", "FISHER_PVALUE_LOG10", bins)

    if "TUMOR_FRACTION" in metric:
        #TUMOR
        bins = np.arange(0, 1, 0.05)

        # graph_tumor_metric(msi_files, "UNIQUE_TUMOR", "UNIQUE_TUMOR_msi", "FISHER_PVALUE_LOG10", bins)
        graph_tumor_metric(mss_files, "UNIQUE_TUMOR", "UNIQUE_TUMOR_mss_more_bunched", "UNIQUE_TUMOR_FRACTION", bins)
        # graph_tumor_metric(gib_files_2, "UNIQUE_TUMOR", "UNIQUE_TUMOR_HG002", "UNIQUE_TUMOR_FRACTION", bins)
        # graph_tumor_metric(gib_files_3, "UNIQUE_TUMOR", "UNIQUE_TUMOR_HG003", "UNIQUE_TUMOR_FRACTION", bins)

    # second_third_normal_allele_counts(gib_files_2, "HG002", old=False)
    # second_third_normal_allele_counts(gib_files_3, "HG003", old=False)

    # mutated_allele_fraction(gib_files_2, "HG002")
    # mutated_allele_fraction(gib_files_3, "HG003")
    #
    # mutated_allele_fraction(mss_files, "mss")
    # mutated_allele_fraction(msi_files, "msi")


if __name__ == '__main__':
    main()
