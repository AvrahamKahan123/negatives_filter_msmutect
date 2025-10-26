import pickle
from typing import List, Tuple

import pandas as pd


def locus_is_heterozygous(reference: int, repeats: List[Tuple[int, int]]) -> bool:
    for rep in repeats:
        for allele in rep:
            if allele!=int(reference) and allele != 0:
                return True
    return False


class PickleReader:
    def __init__(self, fp: str):
        self.pickle_handle = open(fp, 'rb')

    def __next__(self) -> pd.DataFrame:
        try:
            ret = pickle.load(self.pickle_handle)
            return ret
        except EOFError:
            raise StopIteration

    def __iter__(self):
        return self

    def __del__(self):
        self.pickle_handle.close()


def main(fp: str = "loci_population_allele_data.pkl"):
    reader = PickleReader(fp)
    heterozygous_loci = []
    for loci_set in reader:
        print(loci_set.columns)
        print("REFERENCE_REPEATS" in loci_set.columns)
        loci_set["HETEROZYGOUS"] = loci_set.apply(lambda row: locus_is_heterozygous(row["REFERENCE_REPEATS"], [row["ALLELES_REPEATS"]]), axis=1)
        heterozygous_loci.append(loci_set[loci_set['HETEROZYGOUS']])
    all_heterozygous_loci = pd.concat(heterozygous_loci, ignore_index=True)
    all_heterozygous_loci.to_csv("all_heterozygous_loci.csv")


if __name__ == '__main__':
    main()
