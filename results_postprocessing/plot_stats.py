import os.path
from collections import defaultdict
from typing import List, Tuple, Dict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch

from results_postprocessing.cancer_data import results_directory, graphs_directory
from results_postprocessing.enums import COLUMN
from SamplesDB import MSI_CLASSIFICATION


def interleave_list(lst: List[object], interleave_element: str):
    ret = []
    for x in lst:
        ret.append(x)
        ret.append(interleave_element)
    return ret


def split_mutations_df(mutations_df: pd.DataFrame, uuid_column: str, mutations_column: str):

    columns = mutations_df.columns[3:]
    all_tables = []
    for col in columns:
        new_subtable = []
        for _, row in mutations_df.iterrows():
            new_row = {COLUMN.CASE: row[COLUMN.CASE], uuid_column: row[COLUMN.CLASSIFICATION]+"_"+str(col), mutations_column: row[col] }
            new_subtable.append(new_row)
        all_tables.append(new_subtable)
    return pd.concat([pd.DataFrame(subtable) for subtable in all_tables])


# def extract_palette(mutations_df: pd.DataFrame, uuid_column: str, msi_color: str, mss_color: str):
#     pallete = dict()
#     uuids = list(dict(mutations_df[uuid_column].value_counts()).keys())
#     for uuid in uuids:
#         if "MSI" in uuid:
#             pallete[uuid] = msi_color
#         elif "MSS" in uuid:
#             pallete[uuid] = mss_color
#         else:
#             raise RuntimeError(f"Could not classify {uuid} as either mss or msi")
#     return pallete



class UUIDdb:
    def __init__(self,uuid_column = None):
        self.db = defaultdict(dict)
        if uuid_column is not None:
            uuids = list(dict(uuid_column.value_counts()).keys())
            for uuid in uuids:
                self.add_uuid(uuid)

    def add_uuid(self, uuid: str):
        stripped_uuid, ms_status_str = self.split_uuid(uuid)
        self.db[stripped_uuid][ms_status_str] = uuid


    #TODO: why did I strip them?
    def split_uuid(self, uuid: str) -> Tuple[str, str]:
        if "MSI" in uuid:
            stripped_uuid = uuid.replace("MSI", "").replace("_", "")
            category_str = MSI_CLASSIFICATION.MSI
        elif "MSS" in uuid:
            stripped_uuid = uuid.replace("MSS", "").replace("_", "")
            category_str = MSI_CLASSIFICATION.MSS
        elif "NEGATIVE_CONTROL" in uuid:
            stripped_uuid = "GIB"
            category_str = MSI_CLASSIFICATION.NEGATIVE_CONTROL
        else:
            raise RuntimeError(f"Could not classify {uuid} as either mss or msi")

        return stripped_uuid, category_str


    @staticmethod
    def assign_dict_no_overwrite(d: dict, key: object, value: object) -> None:
        if key in d and d[key] is not None:
            raise RuntimeError("Tried to overwrite value")
        else:
            d[key] = value

    def get_palette(self, msi_color: str, mss_color: str) -> Dict[str, str]:
        pallete = dict()
        for category in self.db.keys():
            pallete[self.db[category][MSI_CLASSIFICATION.MSI]] = msi_color
            pallete[self.db[category][MSI_CLASSIFICATION.MSS]] = mss_color
        return pallete

    def uuid_pairs(self, category_order: List[str] = None) -> List[Tuple[str, str]]:
        if category_order is None:
            category_order = list(self.db.keys())
        ret = []
        for category in category_order:
            category = category.replace("_", "")
            mss_uuid = self.db[category][MSI_CLASSIFICATION.MSS]
            msi_uuid = self.db[category][MSI_CLASSIFICATION.MSI]
            ret.append((mss_uuid, msi_uuid))
        return ret


def min_minus_max(mutations_df: pd.DataFrame, uuid_column: str, mutations_column: str, uuid_pairs: List[Tuple[str, str]]) -> List[float]:
    ret = []
    for pair in uuid_pairs:
        if pair[0] =="MSS_0.49":
            croc=1
        mss_part = mutations_df[mutations_df[uuid_column]==pair[0]]
        max_mss = mss_part[mutations_column].max()
        msi_part = mutations_df[mutations_df[uuid_column]==pair[1]]
        min_msi = msi_part[mutations_column].min()
        ret.append(min_msi-max_mss)
    return ret



def create_single_plot(mutations_df: pd.DataFrame, output_dir: str, stat_name: str):
    mutations_df.insert(2, COLUMN.CANCER_TYPE, mutations_df.pop(COLUMN.CANCER_TYPE))
    columns = mutations_df.columns
    num_cols = len(columns) - 3
    cancer_type = mutations_df[COLUMN.CANCER_TYPE].iloc[0]
    uuid_column = "UUID"
    mutations_column = "MUTATIONS"
    mss_color = "skyblue"
    msi_color = "salmon"

    only_has_mss = len(mutations_df["CLASSIFICATION"].value_counts()) == 1
    if only_has_mss: # has only mss
        new_row = pd.DataFrame([{COLUMN.CASE: "FAKE", COLUMN.CANCER_TYPE: cancer_type, COLUMN.CLASSIFICATION: "MSI"} | {col: 0 for col in columns[3:]}])
        mutations_df = pd.concat([mutations_df, new_row], ignore_index=True)

    mutations_df_per_pval = split_mutations_df(mutations_df, uuid_column, mutations_column)
    uuids_manager = UUIDdb(mutations_df_per_pval[uuid_column])
    uuid_pairs = uuids_manager.uuid_pairs(list(mutations_df.columns[3:]))
    palette = uuids_manager.get_palette(msi_color=msi_color, mss_color=mss_color)
    mutations_df_per_pval[mutations_column] = np.log10(np.clip(mutations_df_per_pval[mutations_column], 1, 1e50)) # no clip from above
    fig, ax = plt.subplots(figsize=(num_cols*2, 10))
    separation_metric = min_minus_max(mutations_df_per_pval, uuid_column, mutations_column, uuid_pairs)
    sns.violinplot(data=mutations_df_per_pval, x=uuid_column, y=mutations_column, inner="quartile", ax=ax, palette=palette)
    for i in range(num_cols):
        ax.axvline(x=2*i+1.5, color="black", linewidth=2) # separate between different pvalues
        if not only_has_mss:
            plt.text(x=2*i, y=1.1, s=f"MM: {round(separation_metric[i], 3)}",
                 bbox=dict(facecolor="white", alpha=0.7, edgecolor="black"))

    plt.xlabel(f"{stat_name} THRESHOLD", fontsize=16)
    plt.ylabel("NUMBER OF MUTATIONS (LOG10 SCALE)", fontsize=16)
    ax.set_ylim(1, 7)

    plt.yticks(list(range(1, 8)))

    plt.title(f"{cancer_type}: {stat_name}", fontsize=18)
    plt.tick_params(axis="x", labelsize=14)
    plt.tick_params(axis="y", labelsize=14)

    legend_elements = [
        Patch(facecolor=mss_color, edgecolor="black", label="MSS"),
        Patch(facecolor=msi_color, edgecolor="black", label="MSI")
    ]
    ax.legend(handles=legend_elements, title="Groups")

    xlabels = interleave_list(list(columns[3:]), interleave_element="")
    ax.set_xticklabels(xlabels)
    save_path = f"{cancer_type}_{stat_name}.png"
    print(f"SAVED: {save_path}")
    plt.savefig(os.path.join(output_dir, save_path), dpi=400)
    plt.close(fig)
    # plt.show()


def normal_support_filter(df: pd.DataFrame):
    thresholds = list(range(6, 21))
    df_filtered = df.iloc[:, :3]
    highest_threshold_idx = df.columns.get_loc(str(thresholds[-1]))
    cols = [df.iloc[:, highest_threshold_idx:].apply(sum, axis=1)]
    for current_threshold_idx in range(highest_threshold_idx-1, highest_threshold_idx-len(thresholds), -1): # we already handled last column
        cols.append(df.iloc[:, current_threshold_idx]+cols[-1])
    df_filtered["ORIGINAL"] = df.iloc[:, 3:].apply(sum, axis=1)
    for thresh, col, in zip(thresholds, reversed(cols)):
        df_filtered[str(thresh)] = col
    return df_filtered


def plot_from_csv(csv_fp: str, stat_name: str, plot_negatives: bool = True, output_dir = None, transform_function=None):
    if output_dir is None:
        output_dir = f"{os.path.join(graphs_directory(), stat_name)}"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    results = pd.read_csv(csv_fp)
    if transform_function is not None:
        results = transform_function(results)
    columns = results.columns
    num_cols = len(columns)-3
    plt.figure(figsize=(num_cols*10, 22))
    cancer_types = dict(results[COLUMN.CANCER_TYPE].value_counts())

    if plot_negatives:
        negatives_df = results[results[COLUMN.CANCER_TYPE]==MSI_CLASSIFICATION.NEGATIVE_CONTROL]
    else:
        negatives_df = pd.DataFrame(columns=results.columns)  # empty

    for cancer in cancer_types.keys():
        relevant_df = results[results[COLUMN.CANCER_TYPE] == cancer]
        relevant_df_wnegatives = pd.concat([relevant_df, negatives_df])
        create_single_plot(relevant_df_wnegatives, output_dir, stat_name)


def plot_sample_set_scatter_plot(title: str, save_name: str, set_results: pd.DataFrame, x_column_label: str):
    columns = set_results.columns
    num_cols = len(columns)-3
    plt.figure(figsize=(30, 18))

    for i, col in enumerate(columns[3:]):
        filtered = set_results[set_results["CLASSIFICATION"] == MSI_CLASSIFICATION.MSI]
        msi_column = np.nan_to_num(np.log10(filtered[col]), 1)
        x_msi = np.full(msi_column.shape[0], i)  # now groups are spread along the x-axis

        filtered = set_results[set_results["CLASSIFICATION"] == MSI_CLASSIFICATION.MSS]
        mss_column = np.nan_to_num(np.log10(filtered[col]), 1)  # data along the y-axis
        x_mss = np.full(mss_column.shape[0], i)  # now groups are spread along the x-axis

        filtered = set_results[set_results["CLASSIFICATION"] == MSI_CLASSIFICATION.NEGATIVE_CONTROL]
        negative_column = np.nan_to_num(np.log10(filtered[col]), 1)  # data along the y-axis
        x_negative = np.full(negative_column.shape[0], i)  # now groups are spread along the x-axis

        # three sets of colors per group
        plt.scatter(x_msi, msi_column, c="red", alpha=0.7)
        plt.scatter(x_mss, mss_column, c="blue", alpha=0.7)
        plt.scatter(x_negative, negative_column, c="green", alpha=0.7)

    # Dummy scatters for legend
    plt.scatter([], [], c="red", label="MSI")
    plt.scatter([], [], c="green", label="MSS")
    plt.scatter([], [], c="blue", label="NEGATIVE")

    plt.xticks(range(num_cols), columns[3:])
    plt.xlabel(f"{x_column_label}", fontsize=16)
    plt.ylabel("NUMBER OF MUTATIONS (LOG10 SCALE)", fontsize=16)
    plt.yticks(list(range(1,8)))
    plt.title(title, fontsize=18)
    plt.tick_params(axis="x", labelsize=14)
    plt.tick_params(axis="y", labelsize=14)

    plt.legend(title="Classifications", fontsize=18, title_fontsize=18)
    plt.savefig(save_name, dpi=300)


def normalize(df: pd.DataFrame):
    arr = df.iloc[:, 3:].sum(axis=1).to_numpy()
    return (df.iloc[:, 3:].div(arr, axis=0))


def convert_csv_to_line(csv_fp: str, classification: str) -> List[List[float]]:
    df = pd.read_csv(csv_fp)
    filtered_df = (df[df[COLUMN.CLASSIFICATION]==classification])
    normalized_df = normalize(filtered_df)
    return [list(row.iloc[:]) for _, row in normalized_df.iterrows()]


def convert_csv_to_dist_of_highest_column(csv_fp: str, classification: str, rounding_factor=None) -> List[float]:
    df = pd.read_csv(csv_fp)
    filtered_df = (df[df[COLUMN.CLASSIFICATION] == classification])
    normalized_df = normalize(filtered_df)
    ret = list(normalized_df.iloc[:, -1])
    if rounding_factor is None:
        return ret
    else:
        return [round(x, rounding_factor) for x in ret]


def normalize_histogram(lst: List[float], rounding_factor: int):
    jumps = 10**(-rounding_factor)
    edges = np.arange(0, 1, jumps)
    counts, edges = np.histogram(lst, bins=edges)
    counts = counts/counts.sum()
    return counts, edges


def column_headers(csv_fp: str) -> List[str]:
    df = pd.read_csv(csv_fp)
    return list(df.columns)


def plot_line_graphs(tcga_csv: str, gib_csv: str, title: str, xlabel: str) -> None:
    tcga_mss_lines = convert_csv_to_line(tcga_csv, MSI_CLASSIFICATION.MSS)
    tcga_msi_lines = convert_csv_to_line(tcga_csv, MSI_CLASSIFICATION.MSI)
    gib_lines = convert_csv_to_line(gib_csv, MSI_CLASSIFICATION.NEGATIVE_CONTROL)
    first_col = 3
    headers = column_headers(tcga_csv)
    if "." in headers[first_col+1]:
        transform = float
    else:
        transform = int
    int_columns = [transform(s) for s in headers[first_col:]]
    for dataset in [ (tcga_mss_lines, "blue"), (tcga_msi_lines, "red"),  (gib_lines, "green")]:
    # for dataset in [(tcga_mss_lines, "red")]:

        for sample in dataset[0]:
            plt.plot(int_columns, sample, color=dataset[1])
    # plt.xticks(int_columns[::5])
    plt.ylabel("Mutation Fraction")
    plt.xlabel(xlabel)
    legend_elements = [
        Patch(facecolor="red", label="MSI"),
        Patch(facecolor="blue", label="MSS"),
        Patch(facecolor="green", label="GIB")

    ]

    plt.legend(handles=legend_elements, title="Legend")
    plt.title(title)
    plt.show()


def plot_highest_column_of_line_graph(tcga_csv: str, gib_csv: str, title: str, xlabel: str) -> None:
    tcga_mss_hist, edges = normalize_histogram(convert_csv_to_dist_of_highest_column(tcga_csv, MSI_CLASSIFICATION.MSS, rounding_factor=2), 2)
    tcga_msi_hist, edges = normalize_histogram(convert_csv_to_dist_of_highest_column(tcga_csv, MSI_CLASSIFICATION.MSI, rounding_factor=2), 2)
    gib_hist, edges = normalize_histogram(convert_csv_to_dist_of_highest_column(gib_csv, MSI_CLASSIFICATION.NEGATIVE_CONTROL, rounding_factor=2), 2)
    # plt.bar(bin_edges[:-1], counts, width=bin_widths, align="edge", edgecolor="black")
    plt.bar(edges[:-1], tcga_msi_hist, width=1e-2, align="edge",
            alpha=0.5, color="red", label="MSI")
    plt.bar(edges[:-1], tcga_mss_hist, width=1e-2, align="edge",
            alpha=0.5, color="blue", label="MSS")
    plt.bar(edges[:-1], gib_hist, width=1e-2, align="edge",
            alpha=0.5, color="green", label="GIB")


    plt.xlabel(xlabel)
    plt.ylabel("Percentage of Samples")
    plt.title(title)
    plt.legend()
    plt.show()


def plot_distributions(csv_fp: str, classification: str):
    mutations_df = pd.read_csv(csv_fp)
    mutations_df=mutations_df[mutations_df[COLUMN.CLASSIFICATION]==classification]
    edges = [int(col) for col in mutations_df.columns[3:]]
    for row in list(mutations_df.iterrows())[::10]:
        ns = np.array(list(row[1])[3:])
        plt.plot(edges, ns/ns.sum(),  color="red")
        plt.xlabel("Normal Support")
        plt.ylabel("Mutation Fraction")
        plt.title(classification)
        plt.yticks(np.arange(0, 1.1, 0.1))
    plt.show()


def main():
    plot_from_csv(csv_fp=os.path.join(results_directory(), "full_tcga_NOISELESS_FILTER_NO_UTF.csv"), stat_name="FILTER_COMP")
    # plot_line_graphs("results/full_tcga_locus_length.csv", "results/gib_locus_length.csv", "Locus Length Mutation Fraction")
    # plot_highest_column_of_line_graph("results/full_tcga_purity.csv", "results/gib_purity.csv", "Percentage of Mutations that are Pure", xlabel="Pure Percentage of Mutations")
    # plot_distributions("results/full_tcga_normal_support.csv", classification=MSI_CLASSIFICATION.MSI)
    # plot_line_graphs("results/tcga_test_set_NON_ALLELE_TUMOR_FRACTION.csv",
    #                  "results/gib_NON_ALLELE_TUMOR_FRACTION.csv",
    #                  "NON_ALLELE_TUMOR_FRACTION", "Non-Allele Tumor Fraction")
    # plot_distributions("results/full_tcga_normal_support.csv", classification=MSI_CLASSIFICATION.MSS)


if __name__ == '__main__':
    main()