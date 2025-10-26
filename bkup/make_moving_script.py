import os, shutil
from typing import List, Set, Tuple


def classifications(fp: str = "C:/Users/avrah/MaruvkaLab/Texas_samples_organization/tcga_msi_classification_FULL_11_09_25.txt") -> Tuple[Set[str], Set[str]]:
    # returns the MSI, MSS samples lists as sets
    with open(fp, 'r') as croc:
        all_lines = croc.readlines()
    mss_samps = set()
    msi_samps = set()
    for l in all_lines:
        split_line = l.split("\t")
        classification = split_line[1].rstrip()
        if classification == "MSS":
            mss_samps.add(f"{split_line[0]}")
        elif classification == "MSI":
            msi_samps.add(f"{split_line[0]}")
    return msi_samps, mss_samps


def find_first_period_or_underscore(s: str) -> int:
    """
    Find the index of the first period (.) or underscore (_) in a string.
    """
    # Find index of '.' and '_'
    period_index = s.find('.')
    underscore_index = s.find('_')

    # If one of them is -1, return the other
    if period_index == -1:
        return underscore_index
    if underscore_index == -1:
        return period_index

    # Return the earlier occurrence
    return min(period_index, underscore_index)


def sample_name(filename: str) -> str:
    basename = os.path.basename(filename)
    return basename[:find_first_period_or_underscore(basename)]


def main():
    msi, mss = classifications()
    src_dir = "/data/MSMuTect_called_mut_filt_fixed/"
    all_files = [os.path.join(src_dir, file) for file in os.listdir(src_dir)]
    for f in all_files:
        sample = sample_name(f)
        if sample in msi:
            dest_path = "/data/msi"
        elif sample in mss:
            dest_path = "/data/mss"
        else:
            raise RuntimeError(f"Unknown Sample: {f}")

        shutil.move(f, dest_path)


if __name__ == '__main__':
    main()