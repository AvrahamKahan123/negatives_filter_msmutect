import glob, os, time

import numpy as np

LOCI_LINE_COUNT = 27705407


def all_npy_files(src_dir: str):
    return glob.glob(os.path.join(src_dir, "*.npy"))


def main():
    npy_files = all_npy_files('/storage/bfe_maruvka/avrahamk/Negatives/population_level_analysis')
    final_file_fp = "/storage/bfe_maruvka/avrahamk/Negatives/final_pop.npy"
    combined_data = np.zeros((LOCI_LINE_COUNT, 41), dtype=np.int32)
    for f in npy_files:
        st = time.time()
        current_data = np.load(f)
        combined_data = combined_data + current_data
        e = time.time()
        print(f"{f}: {e-st}")
    np.save(final_file_fp, combined_data)
