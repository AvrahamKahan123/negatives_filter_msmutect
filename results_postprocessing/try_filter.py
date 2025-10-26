import pandas as pd
from results_postprocessing.analyze_gib import filter_results
from utils import filtered_mutation_files_list


class SampleSet:
    def __init__(self, name: str, origin_dir: str, evaluate: bool):
        self.name = name
        self.file_list = filtered_mutation_files_list(origin_dir)
        self.evaluate = evaluate


def main(mss: bool, msi: bool, gib: bool):
    msi_sample_set = SampleSet("MSI", "../data/msi_example_files", msi)
    mss_sample_set = SampleSet("MSS", "../data/mss_example_files", mss)
    gib_sample_set = SampleSet("GIB", "../data/gib_files_filtered", gib)
    all_dicts = []
    for sample_set in [gib_sample_set, mss_sample_set, msi_sample_set]:
        if sample_set.evaluate:
            print(f"*******************************{sample_set.name}*******************************")
            category = "GIB"
            for file in sample_set.file_list:
                old_muts, post_filt_muts = filter_results(file)
                all_dicts.append(
                    {"FILE": file, "Category": category, "Old_Mutations": old_muts, "New_Mutations": post_filt_muts})
        combined_df = pd.DataFrame(all_dicts)
        combined_df.to_csv("filter_test_3.csv")
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


if __name__ == '__main__':
    main(gib=True, mss=False, msi=False)
