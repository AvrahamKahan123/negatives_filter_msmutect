import operator, os, time, sys
from collections import defaultdict
from dataclasses import dataclass

from typing import Dict, List, Union, Set, Tuple, Callable
import multiprocessing as mp
import numpy as np
import pandas as pd

from cancer_data import load_all_samples_gib, load_all_samples, load_test_samples, data_directory, results_directory
from create_purity_db import DB_connection, Locus, connect_to_purity_db
from enums import COLUMN, MSI_CLASSIFICATION
from SamplesDB import SamplesDB, Sample, SamplesSet
from results_postprocessing.NoisyLocusDB import NoisyLocusDB


def column_headers(prefix: str, thresholds: List[Union[float, int]]) -> List[str]:
        return [f"{prefix}{thresh}" for thresh in thresholds]


def fisher_thresholds(sample: Sample) -> Dict[str, object]:
    log_thresholds = np.flip(np.arange(-5, -1.5, 0.25), axis=0)
    fisher_column_headers = column_headers("P=", log_thresholds)
    # fisher_headers = ["CASE", "CLASSIFICATION", "ORIGINAL"]
    mutations_df = pd.read_csv(sample.mutations_tsv_fp, delimiter="\t")
    thresholds = list(10 ** np.array(log_thresholds))
    mutations_per_threshold = ({"CASE": sample.sample_name(), "CLASSIFICATION": sample.classification, "ORIGINAL": len(mutations_df)} |
                               {thresh_str: (mutations_df[COLUMN.FISHER_TEST_P_VALUE] <= thresh).sum() for thresh_str, thresh in
                                                     zip(fisher_column_headers, thresholds)})
    return mutations_per_threshold


def second_third_normal_allele_counts(fps, title, old=True):
    print(f"*************\n{title}")
    for fp in fps:
        print(os.path.basename(fp))
        df = pd.read_csv(fp, delimiter="\t")
        num_muts = len(df)
        print(f"{len(df[df['NORMAL_MOTIF_REPEATS_2']!=0])/num_muts}, {len(df[df['NORMAL_MOTIF_REPEATS_3']!=0])/num_muts}")
        # print(f"NUMBER OF SECOND ALLELES: {len(df[df['NORMAL_MOTIF_REPEATS_2']!=0])}")
        # print(f"NUMBER OF THIRD ALLELES: {len(df[df['NORMAL_MOTIF_REPEATS_3']!=0])}")


def calc_unique_allele_fraction(normal_motifs: Set[int], tumor_motifs: Dict[int, int]) -> float:
    unique = 0
    total = sum(tumor_motifs.values())
    for motif in tumor_motifs.keys():
        if motif not in normal_motifs:
            unique+=tumor_motifs[motif]
    return unique/total


def add_mutated_allele_fraction_column(df: pd.DataFrame) -> List[pd.DataFrame]:
    df["UNIQUE_TUMOR_FRACTION"] = df.apply(lambda row: calc_unique_allele_fraction(
        set([row[f"NORMAL_MOTIF_REPEATS_{i}"] for i in range(1, 6)]),
        {row[f"TUMOR_MOTIF_REPEATS_{i}"]: row[f"TUMOR_SUPPORTING_READS_{i}"] for i in range(1, 6)}),
        axis=1)


def unique_tumor_fraction(sample: Sample) -> Dict[str, object]:
    st = time.time()
    thresholds = [round(x, 2) for x in list(np.arange(0.01, 0.5, 0.03))]
    utf_headers = column_headers("", thresholds)
    mutations_df = pd.read_csv(sample.mutations_tsv_fp, delimiter="\t")

    add_mutated_allele_fraction_column(mutations_df)
    mutations_per_threshold = (
                {"CASE": sample.sample_name(), "CLASSIFICATION": sample.classification, "ORIGINAL": len(mutations_df)} |
                {thresh_str: (mutations_df["UNIQUE_TUMOR_FRACTION"] >= thresh).sum() for thresh_str, thresh in
                 zip(utf_headers, thresholds)})
    e=time.time()
    # mutations_df.to_csv(f"data/gib_files_wtumorfraction/{sample.sample_name()}.csv")
    print(f"Num Rows: {len(mutations_df)}: {e-st}")
    return mutations_per_threshold


def ks_thresholds(sample: Sample) -> Dict[str, object]:
    log_thresholds = np.flip(np.arange(-7, 0, 0.25), axis=0)
    fisher_column_headers = column_headers("P=", log_thresholds)
    mutations_df = pd.read_csv(sample.mutations_tsv_fp, delimiter="\t")
    thresholds = list(10 ** np.array(log_thresholds))
    mutations_per_threshold = ({"CASE": sample.sample_name(), "CLASSIFICATION": sample.classification, "ORIGINAL": len(mutations_df)} |
                               {thresh_str: (mutations_df["KS_TEST_PVALUE"] <= thresh).sum() for thresh_str, thresh in
                                                     zip(fisher_column_headers, thresholds)})
    return mutations_per_threshold


def separate_keys_into_regular_and_oversize(d: dict, max_val: int, from_above: bool) -> Tuple[list, list]:
    regular = []
    overlarge = []
    if from_above:
        comparison = operator.ge
    else:
        comparison = operator.le

    for key in d.keys():
        if comparison(key, max_val):
            overlarge.append(key)
        else:
            regular.append(key)
    return regular, overlarge


def bound_values(d: dict, max_val: int, from_above=True) -> defaultdict:
    # reduces all values greater than max val to max value
    ret = defaultdict(int)
    regular, overlarge_keys = separate_keys_into_regular_and_oversize(d, max_val, from_above)
    for key in regular:
        ret[key] = d[key]
    for key in overlarge_keys:
        ret[max_val]+=d[key]
    return ret


def motif_length(sample: Sample) -> Dict[str, object]:
    min_size = 1
    max_size = 15
    lengths = list(range(min_size, max_size))
    length_column_headers = column_headers("", lengths)
    mutations_df = pd.read_csv(sample.mutations_tsv_fp, delimiter="\t")
    mutations_df["motif_length"] = mutations_df["PATTERN"].apply(len)
    value_counts = dict(mutations_df["motif_length"].value_counts())
    value_counts_bounded = bound_values(value_counts, max_size)
    mutations_per_threshold = ({"CASE": sample.sample_name(), "CLASSIFICATION": sample.classification} |
                               {str(i): value_counts_bounded[i] for i in range(min_size, max_size + 1)})
    return mutations_per_threshold


def allele_column_name(df: pd.DataFrame) -> str:
    # this is to cover for a bad typo error in MSMuTect <= 4.1
    cols = df.columns
    if "TUMOR_ALLELES_2" in cols:
        return "ALLELES"
    else:
        return "ALLELE"


def add_non_allele_tumor_fraction_column(df: pd.DataFrame):
    df[COLUMN.NON_ALLELE_TUMOR_FRACTION] = df.apply(lambda row: calc_unique_allele_fraction(
        set([row[f"TUMOR_ALLELE_1"]]+[row[f"TUMOR_{allele_column_name(df)}_{i}"] for i in range(2, 5)]),
        {row[f"TUMOR_MOTIF_REPEATS_{i}"]: row[f"TUMOR_SUPPORTING_READS_{i}"] for i in range(1, 6)}),
                                           axis=1)


def distribution_histogram(df: pd.DataFrame, column: str, rounding_factor: int = None):
    if rounding_factor is not None:
        df[column] = round(df[column], rounding_factor)
    counts = dict(df[column].value_counts())
    bins = defaultdict(int)
    for val, count in counts.items():
        bins[val] = count
    return bins


def non_allele_tumor_fraction(sample) -> Dict[str, object]:
    thresholds = [round(x, 2) for x in list(np.arange(0.01, 1.01, 0.01))]
    utf_headers = column_headers("", thresholds)
    mutations_df = pd.read_csv(sample.mutations_tsv_fp, delimiter="\t")
    add_non_allele_tumor_fraction_column(mutations_df)
    bins = distribution_histogram(mutations_df, COLUMN.NON_ALLELE_TUMOR_FRACTION, 2)
    # add_mutated_allele_fraction_column(mutations_df)
    mutations_per_threshold = (
                {"CASE": sample.sample_name(), "CLASSIFICATION": sample.classification} |
                {thresh_str: bins[thresh] for thresh_str, thresh in zip(utf_headers, thresholds)})
    return mutations_per_threshold


def normal_support(sample: Sample) -> Dict[str, object]:
    min_support = 5
    max_support = 75
    lengths = list(range(min_support, max_support))
    mutations_df = pd.read_csv(sample.mutations_tsv_fp, delimiter="\t")
    mutations_df["NORMAL_S"] = (mutations_df["NORMAL_SUPPORTING_READS_1"]+mutations_df["NORMAL_SUPPORTING_READS_2"]+
                                mutations_df["NORMAL_SUPPORTING_READS_3"]+mutations_df["NORMAL_SUPPORTING_READS_4"]+
                                mutations_df["NORMAL_SUPPORTING_READS_5"])
    value_counts = dict(mutations_df["NORMAL_S"].value_counts())
    value_counts_bounded = bound_values(value_counts, max_support)
    mutations_per_threshold = ({"CASE": sample.sample_name(), "CLASSIFICATION": sample.classification} |
                               {str(i): value_counts_bounded[i] for i in range(min_support, max_support + 1)})
    return mutations_per_threshold


def locus_length(sample: Sample) -> Dict[str, object]:
    min_size = 6
    max_size = 100
    lengths = list(range(min_size, max_size))
    length_column_headers = column_headers("", lengths)
    mutations_df = pd.read_csv(sample.mutations_tsv_fp, delimiter="\t")
    mutations_df["locus_length"] = mutations_df["REFERENCE_SEQUENCE"].apply(len)
    value_counts = dict(mutations_df["locus_length"].value_counts())
    value_counts_bounded = bound_values(value_counts, max_size)
    mutations_per_threshold = ({"CASE": sample.sample_name(), "CLASSIFICATION": sample.classification} |
                               {str(i): value_counts_bounded[i] for i in range(min_size, max_size+1)})
    return mutations_per_threshold


def purity_func(row, db_connection: DB_connection):
    return db_connection.query_locus_purity(Locus(f'chr{row["CHROMOSOME"]}', row["START"], row["END"], row["PATTERN"], None))


def purity_filter(sample: Sample, db_connection: DB_connection) -> Dict[str, object]:
    min_val = 85
    max_val = 100
    # lengths = list(range(min_val, max_val+1))
    mutations_df = pd.read_csv(sample.mutations_tsv_fp, delimiter="\t")

    mutations_df["purity"] = mutations_df.apply(purity_func, axis=1, db_connection=db_connection)

    # mutations_df["purity"] = mutations_df["REFERENCE_SEQUENCE"].apply(locus_purity)
    value_counts = dict(mutations_df["purity"].value_counts())
    value_counts_bounded = bound_values(value_counts, min_val, from_above=False)
    mutations_per_threshold = ({"CASE": sample.sample_name(), "CLASSIFICATION": sample.classification} |
                               {str(i): value_counts_bounded[i] for i in range(min_val, max_val+1)})
    return mutations_per_threshold


def split_keys_by_threshold(keys, thresholds, comparator) -> List[List[int]]:
    return [[k for k in keys if comparator(k, thresh)] for thresh in thresholds]


def evaluate_passes_each_threshold(df: pd.DataFrame, column: str, thresholds: List[int], comparator):
    # should only be used for int thresholds
    val_counts: Dict[int, int] = dict(df[column].value_counts())
    relevant_keys = split_keys_by_threshold(val_counts.keys(), thresholds, comparator)
    num_mutations = []
    for key_set in relevant_keys:
        current_num_muts = 0
        for key in key_set:
            current_num_muts+=val_counts[key]
        num_mutations.append(current_num_muts)
    return num_mutations


def avraham_filter(mutations_df: pd.DataFrame, normal_support_threshold: int = 8, UTF: float = 0.1,
                   third_and_later_fraction: float = 0.1, fisher_threshold = 1):
    mutations_df["NORMAL_S"] = (mutations_df["NORMAL_SUPPORTING_READS_1"] +
                                mutations_df["NORMAL_SUPPORTING_READS_2"] +
                                mutations_df["NORMAL_SUPPORTING_READS_3"] +
                                mutations_df["NORMAL_SUPPORTING_READS_4"])
    mutations_df["TUMOR_S"] = (mutations_df["TUMOR_SUPPORTING_READS_1"] +
                               mutations_df["TUMOR_SUPPORTING_READS_2"] +
                               mutations_df["TUMOR_SUPPORTING_READS_3"] +
                               mutations_df["TUMOR_SUPPORTING_READS_4"])

    mutations_df["THIRD_AND_LATER_TUMOR_PROPORTION"] = (mutations_df["TUMOR_SUPPORTING_READS_3"] + mutations_df[
        "TUMOR_SUPPORTING_READS_4"]) / mutations_df["TUMOR_S"]
    mutations_df["THIRD_AND_LATER_NORMAL_PROPORTION"] = (mutations_df["NORMAL_SUPPORTING_READS_3"] + mutations_df[
        "NORMAL_SUPPORTING_READS_4"]) / mutations_df["NORMAL_S"]
    mutations_df["UNIQUE_TUMOR_FRACTION"] = mutations_df.apply(lambda row: calc_unique_allele_fraction(
        set([row[f"NORMAL_MOTIF_REPEATS_{i}"] for i in range(1, 6)]),
        {row[f"TUMOR_MOTIF_REPEATS_{i}"]: row[f"TUMOR_SUPPORTING_READS_{i}"] for i in range(1, 6)}),
                                                               axis=1)

    mutations_df_filtered = mutations_df[
        (mutations_df["NORMAL_S"] >= normal_support_threshold) &
        (mutations_df["UNIQUE_TUMOR_FRACTION"] > UTF) &
        (mutations_df["THIRD_AND_LATER_NORMAL_PROPORTION"] < third_and_later_fraction) &
        (mutations_df[COLUMN.FISHER_TEST_P_VALUE] < fisher_threshold)
        ]
    return mutations_df_filtered


def noiseless_filter_without_utf(sample: Sample) -> Dict[str, object]:
    thresholds = list(range(8, 16))
    mutations_df = pd.read_csv(sample.mutations_tsv_fp, delimiter="\t")
    mutations_df_filtered = avraham_filter(mutations_df, UTF=0.0, fisher_threshold=0.01)
    num_mutations_per_threshold = evaluate_passes_each_threshold(mutations_df_filtered, "NORMAL_S", thresholds,
                                                                 operator.ge)
    mutations_per_threshold = (
                {"CASE": sample.sample_name(), "CLASSIFICATION": sample.classification, "ORIGINAL": len(mutations_df)} |
                {str(thresh): muts for thresh, muts in zip(thresholds, num_mutations_per_threshold)})
    return mutations_per_threshold


def noiseless_filter(sample: Sample) -> Dict[str, object]:
    thresholds = list(range(8, 16))
    mutations_df = pd.read_csv(sample.mutations_tsv_fp, delimiter="\t")
    mutations_df_filtered = avraham_filter(mutations_df)
    num_mutations_per_threshold = evaluate_passes_each_threshold(mutations_df_filtered, "NORMAL_S", thresholds,
                                                                 operator.ge)

    mutations_per_threshold = (
                {"CASE": sample.sample_name(), "CLASSIFICATION": sample.classification, "ORIGINAL": len(mutations_df)} |
                {str(thresh): muts for thresh, muts in zip(thresholds, num_mutations_per_threshold)})
    return mutations_per_threshold


def normal_threshold_per_locus(mutations_df: pd.DataFrame, noisy_db: NoisyLocusDB, noiseless_threshold: int, noisy_threshold: int) -> np.array:
    # returns the threshold for each locus based on population level data
    # noisy_data = np.load(os.path.join(data_directory(), "noisy_locus.npy"))
    # return np.where(noisy_data, noisy_threshold, noiseless_threshold)
    thresholds = np.zeros(mutations_df.shape[0], dtype=np.int32)
    for row_num, row in mutations_df.iterrows():
        if noisy_db.locus_present_in_db(row[COLUMN.CHROMOSOME], row[COLUMN.START],
                                        row[COLUMN.STOP], row[COLUMN.PATTERN]):
            current_threshold = noisy_threshold
        else:
            current_threshold = noiseless_threshold
        thresholds[row_num] = current_threshold
    return thresholds

def yossi_filter(mutations_df: pd.DataFrame, noisy_db: NoisyLocusDB) -> pd.DataFrame:
    normal_support_threshold_noiseless = 8
    normal_support_threshold_noisy = 15
    normal_third_motif_plus_threshold = 0.1
    fisher_threshold = 0.01
    # UTF = 0
    mutations_df["NORMAL_S"] = (mutations_df["NORMAL_SUPPORTING_READS_1"] +
                                mutations_df["NORMAL_SUPPORTING_READS_2"] +
                                mutations_df["NORMAL_SUPPORTING_READS_3"] +
                                mutations_df["NORMAL_SUPPORTING_READS_4"] +
                                mutations_df["NORMAL_SUPPORTING_READS_5"])

    mutations_df["THIRD_AND_LATER_NORMAL_PROPORTION"] = (mutations_df["NORMAL_SUPPORTING_READS_3"] + mutations_df[
        "NORMAL_SUPPORTING_READS_4"] + mutations_df["NORMAL_SUPPORTING_READS_5"]) / mutations_df["NORMAL_S"]

    threshold_per_locus = normal_threshold_per_locus(mutations_df, noisy_db, normal_support_threshold_noiseless, normal_support_threshold_noisy)
    print((threshold_per_locus==normal_support_threshold_noisy).mean())
    mutations_df_filtered = mutations_df[
        (mutations_df["NORMAL_S"] > threshold_per_locus) &
        (mutations_df[COLUMN.FISHER_TEST_P_VALUE] <= fisher_threshold) &
        (mutations_df["THIRD_AND_LATER_NORMAL_PROPORTION"] < normal_third_motif_plus_threshold)
        ]

    return mutations_df_filtered


def compare_filters(sample: Sample, noisy_db: NoisyLocusDB) -> Dict[str, object]:
    print(sample.classification)
    mutations_df = pd.read_csv(sample.mutations_tsv_fp, delimiter="\t")
    if sample.classification==MSI_CLASSIFICATION.MSI:
        croc=1
    yossi_filtered = yossi_filter(mutations_df.copy(), noisy_db)
    avraham_filtered = avraham_filter(mutations_df.copy())
    mutations_per_threshold = ({"CASE": sample.sample_name(), "CLASSIFICATION": sample.classification, "ORIGINAL": len(mutations_df)} |
                               {"AVRAHAM_FILTER": len(avraham_filtered), "YOSSI_FILTER": len(yossi_filtered)})
    return mutations_per_threshold


def split_into_work_chunks(samples: SamplesSet, minimum_chunk_size: int = 3e7) -> List[SamplesSet]:
    current_chunk_size = 0
    chunks = []
    current_chunk = []
    for samp in samples.mss_patients+samples.msi_patients+samples.negative_patients:
        current_chunk_size+=os.path.getsize(samp.mutations_tsv_fp)
        current_chunk.append(samp)
        if current_chunk_size > minimum_chunk_size:
            chunk_sample_set = SamplesSet(samples.cancer_type)
            chunk_sample_set.add_multiple_samples(current_chunk)
            chunks.append(chunk_sample_set)
            current_chunk = []
            current_chunk_size = 0
    if len(current_chunk)>0:
        new_sample_set = SamplesSet(samples.cancer_type)
        new_sample_set.add_multiple_samples(current_chunk)
        chunks.append(new_sample_set)
    return chunks


def tabulate_statistic(samples: SamplesSet, func: Callable, args: list) -> pd.DataFrame:
    # function analyzes every sample in a given way. It must do the following: Accept a sample as input and output a dict with at least two columns titled "Case" and "Classification"
    # returns a dataframe with column Case: Sample_name, Classification: classification, and then columns holding the number of mutations under different parameters
    sample_results = []
    for sample in (samples.mss_patients+samples.msi_patients+samples.negative_patients):
        st = time.time()
        statistic_row = func(sample, *args)
        sample_results.append(statistic_row)
        e = time.time()
        print(f"time for sample: {e - st}")
    ret = pd.DataFrame(sample_results)
    print(f"{samples.cancer_type}: {ret.shape}")
    ret.insert(1, COLUMN.CANCER_TYPE, samples.cancer_type)

    return ret


def tabulate_statistic_parallel(samples: SamplesSet, func: Callable, args: list, num_cpus: int = 8) -> pd.DataFrame:
    work_chunks = split_into_work_chunks(samples)
    with mp.Pool(processes=num_cpus) as pool:
        results = pool.starmap(tabulate_statistic, [(chunk, func, args) for chunk in work_chunks])
    return pd.concat(results, ignore_index=True)


class SampleChoice:
    def __init__(self, TCGA: bool, test_dataset_only: bool = False, mss_only: bool = False):
        # True=TCGA. False=GIB
        self.TCGA = TCGA
        self.mss_only = mss_only
        self.test_dataset_only = test_dataset_only

        if self.GIB and self.test_dataset_only:
            raise RuntimeError("Illogical state. GIB has no test dataset")

    @property
    def GIB(self):
        return not self.TCGA

@dataclass
class Test:
    func: Callable
    test_name: str
    args: List[object]


def lookup_test(test_name: str) -> Test:
    tests: Dict[str, Test] = {
        "FISHER": Test(fisher_thresholds, "fisher", []),
        "KS": Test(ks_thresholds, "ks", []),
        "UNIQUE_TUMOR_FRACTION": Test(unique_tumor_fraction, "Unique_Tumor_Fraction", []),
        "LOCUS_LENGTH": Test(locus_length, "locus_length", []),
        "MOTIF_LENGTH": Test(motif_length, "motif_length", []),
        "NORMAL_SUPPORT": Test(normal_support, "normal_support", []),
        "PURITY": Test(purity_filter, "purity", [connect_to_purity_db()]),
        "NOISELESS_FILTER": Test(noiseless_filter, "noiseless_filter", []),
        "NON_ALLELE_TUMOR_FRACTION": Test(non_allele_tumor_fraction, "non_allele_tumor_fraction", []),
        "NOISELESS_FILTER_NO_UTF": Test(noiseless_filter_without_utf, "NOISELESS_FILTER_NO_UTF", [])
        # "COMPARE_FILTERS": Test(compare_filters, "compare_filters", [NoisyLocusDB()]) # not ideal; loads even when going to be used!
    }
    if test_name not in tests:
        raise RuntimeError(f"Unrecognized test: {test_name}\nRecognized tests are {tests.keys()}")
    else:
        return tests[test_name]


def main(test_name: str, sample_choice: SampleChoice, parallel: bool = True):

    test = lookup_test(test_name)

    if sample_choice.GIB:
        all_samples = load_all_samples_gib()
        prefix = "gib"
    else:
        if sample_choice.test_dataset_only:
            all_samples = load_test_samples(mss_only=sample_choice.mss_only)
            prefix = "tcga_test_set"
        else:
            all_samples = load_all_samples(mss_only=sample_choice.mss_only)
            prefix = "full_tcga"

    results = []
    for sample_set in all_samples.get_all_patients():
        if parallel:
            result = tabulate_statistic_parallel(sample_set, test.func, test.args) # parallelization is on this level because most of the runtime is for specific sample sets
        else:
            result = tabulate_statistic(sample_set, test.func, test.args)

        results.append(result)

    all_stats = pd.concat(list(results), ignore_index=True)
    all_stats.to_csv(f"{os.path.join(results_directory(), prefix)}_{test_name}.csv", index=False)


if __name__ == '__main__':
    st = time.time()
    main(test_name="NOISELESS_FILTER_NO_UTF", sample_choice=SampleChoice(TCGA=True, test_dataset_only=False), parallel=True)
    e = time.time()
    print(e-st)
