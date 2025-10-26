from typing import Tuple

import pandas as pd
import os
from bkup import make_moving_script
from results_postprocessing.analyze_full_tcga_mut_set import calc_unique_allele_fraction
from utils import filtered_mutation_files_list


def gib_sample_name(filename: str) -> str:
    basename = os.path.basename(filename)
    return basename[:basename.find(".")]


def sample_name(fp: str):
    if "TCGA" in fp.upper():
        return make_moving_script.sample_name(fp)
    else:
        return gib_sample_name(fp)


def filter_results(fp: str, save: bool = False) -> Tuple[int, int]:
    mutations_df = pd.read_csv(fp, delimiter="\t")
    original_num_mutations = len(mutations_df)
    mutations_df["NORMAL_S"] =            (mutations_df["NORMAL_SUPPORTING_READS_1"] +
                                          mutations_df["NORMAL_SUPPORTING_READS_2"] +
                                          mutations_df["NORMAL_SUPPORTING_READS_3"] +
                                          mutations_df["NORMAL_SUPPORTING_READS_4"])
    mutations_df["TUMOR_S"] = (mutations_df["TUMOR_SUPPORTING_READS_1"] +
                                mutations_df["TUMOR_SUPPORTING_READS_2"] +
                                mutations_df["TUMOR_SUPPORTING_READS_3"] +
                                mutations_df["TUMOR_SUPPORTING_READS_4"])

    mutations_df["THIRD_AND_LATER_TUMOR_PROPORTION"] =  (mutations_df["TUMOR_SUPPORTING_READS_3"] + mutations_df["TUMOR_SUPPORTING_READS_4"])/mutations_df["TUMOR_S"]
    mutations_df["THIRD_AND_LATER_NORMAL_PROPORTION"] = (mutations_df["NORMAL_SUPPORTING_READS_3"] + mutations_df["NORMAL_SUPPORTING_READS_4"])/mutations_df["NORMAL_S"]
    mutations_df["UNIQUE_TUMOR_FRACTION"] = mutations_df.apply(lambda row: calc_unique_allele_fraction(
        set([row[f"NORMAL_MOTIF_REPEATS_{i}"] for i in range(1, 6)]),
        {row[f"TUMOR_MOTIF_REPEATS_{i}"]: row[f"TUMOR_SUPPORTING_READS_{i}"] for i in range(1, 6)}),
                                             axis=1)
    # mutations_df_filtered = mutations_df[((mutations_df["NORMAL_S"] >= 12) | (mutations_df["THIRD_AND_LATER_TUMOR_PROPORTION"]>0.1)) &
    #                                       (mutations_df["UNIQUE_TUMOR_FRACTION"] > 0.1) &
    #                                       (mutations_df["THIRD_AND_LATER_NORMAL_PROPORTION"] < 0.1)
    #                                       ]
    mutations_df_filtered = mutations_df[
        (mutations_df["NORMAL_S"] >= 13) &
        (mutations_df["UNIQUE_TUMOR_FRACTION"] > 0.1) &
        (mutations_df["THIRD_AND_LATER_NORMAL_PROPORTION"] < 0.1)
        ]

    conditions = ["NORMAL SUPPORT", "UNIQUE TUMOR FRACTION", "THIRD_AND_LATER"]
    filters = [(mutations_df["NORMAL_S"] >= 13),
    (mutations_df["UNIQUE_TUMOR_FRACTION"] > 0.1),
    (mutations_df["THIRD_AND_LATER_NORMAL_PROPORTION"] < 0.1)]
    filter_sums = [f.sum()/len(f) for f in filters]
    for cond, f in zip(conditions, filter_sums):
        print(f"{cond}: {f}")

    case_name = sample_name(fp)
    print(f"Case Name: {case_name}\nOriginal Num Mutations: {original_num_mutations}\nFiltered Num Mutations: {len(mutations_df_filtered)}\nPERCENTAGE LEFT: {len(mutations_df_filtered)/original_num_mutations}*****************************")
    if save:
        mutations_df_filtered.to_csv(f"STRICT_{case_name}.csv")
    return original_num_mutations, len(mutations_df_filtered)



if __name__ == '__main__':
    # msi_example_files = filtered_mutation_files_list("data/msi_example_files")
    # mss_example_files = filtered_mutation_files_list("data/mss_example_files")
    gib_files = filtered_mutation_files_list("data/gib_files_filtered")
    #
    all_dicts = []
    print("*******************************GIB*******************************")
    category = "GIB"
    for file in gib_files:
        old_muts, post_filt_muts = filter_results(file)
        all_dicts.append({"FILE": file, "Category": category, "Old_Mutations": old_muts, "New_Mutations": post_filt_muts})
    # print("*******************************MSS*******************************")
    # category = "MSS"
    # for file in mss_example_files:
    #     old_muts, post_filt_muts = filter_results(file)
    #     all_dicts.append({"FILE": file, "Category": category, "Old_Mutations": old_muts, "New_Mutations": post_filt_muts})
    # print("*******************************MSI*******************************")
    # category = "MSI"
    # for file in msi_example_files:
    #     old_muts, post_filt_muts = filter_results(file)
    #     all_dicts.append({"FILE": file, "Category": category, "Old_Mutations": old_muts, "New_Mutations": post_filt_muts})
    # combined_df = pd.DataFrame(all_dicts)
    # combined_df.to_csv("filter_test_2.csv")
    # main(load_gib=True, param="PURITY")
    # tmp = pd.read_csv("data/gib_files_filtered/hg001_normal2_tumor1.filtered.full.mut.tsv", delimiter="\t")
    # tmp = pd.read_csv("data/msi/TCGA-AX-A0IZ.called.filt.mut.tsv.gz", delimiter="\t")
    # print(len(tmp))
    # tmp["NORMAL_S"] = tmp["NORMAL_SUPPORTING_READS_1"]+tmp["NORMAL_SUPPORTING_READS_2"]+tmp["NORMAL_SUPPORTING_READS_3"]+tmp["NORMAL_SUPPORTING_READS_4"]
    # tmp["UNIQUE_TUMOR_FRACTION"] = tmp.apply(lambda row: calc_unique_allele_fraction(
    #     set([row[f"NORMAL_MOTIF_REPEATS_{i}"] for i in range(1, 6)]),
    #     {row[f"TUMOR_MOTIF_REPEATS_{i}"]: row[f"TUMOR_SUPPORTING_READS_{i}"] for i in range(1, 6)}),
    #                                        axis=1)
    # tmp = tmp[((tmp["NORMAL_S"]>12) & (tmp["UNIQUE_TUMOR_FRACTION"]>0.1) & (tmp["NORMAL_MOTIF_REPEATS_3"]==0))]
    # tmp.to_csv("high_support.csv", index=False)
    # print(len(tmp))
    # main(load_gib=True, param="NOISELESS_FILTER")
