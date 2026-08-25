import pandas as pd

from results_postprocessing.analyze_full_tcga_mut_set import avraham_filter
from results_postprocessing.cancer_data import load_all_samples_gib


def main():
    length_counts = {x: 0 for x in range(1, 30)}
    all_samples = load_all_samples_gib()
    for patient in all_samples.get_all_patients()[0].negative_patients:
        print(patient)
        df = pd.read_csv(patient.mutations_tsv_fp, sep="\t")
        filtered = avraham_filter(df, normal_support_threshold=10, UTF=0.0, third_and_later_fraction=0.1,fisher_threshold = 0.01).copy()
        filtered["pattern_length"] = filtered["PATTERN"].str.len()

        for length, count in filtered["pattern_length"].value_counts().items():
            length_counts[length] = length_counts.get(length, 0) + count
    print(length_counts)
    print(sum(length_counts.values()))

if __name__ == "__main__":
    main()