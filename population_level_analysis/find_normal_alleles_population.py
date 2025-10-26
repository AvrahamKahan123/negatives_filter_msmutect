import glob, zipfile, os, sys, time, gzip
from typing import List
import numpy as np
import pandas as pd

LOCI_LINE_COUNT = 27705407


def safe_convert_to_int(s: str):
    if s=="NA":
        return 0
    else:
        return int(s)


def create_matrix_for_tsv(msmutect_fp: str, matrix: np.ndarray):
    # changes the matrix inplace
    st = time.time()
    with open(msmutect_fp, 'r') as msmutect_output:
    # with gzip.open(msmutect_fp, "rt") as msmutect_output:  # "rt" = read text
        header_line = msmutect_output.readline()
        allele_1_col = 20
        allele_2_col = allele_1_col+1
        allele_3_col = allele_2_col+1
        for line_num, line in enumerate(msmutect_output):
            split_line = line.split("\t")
            normal_allele_1 = safe_convert_to_int(split_line[allele_1_col])
            normal_allele_2 = safe_convert_to_int(split_line[allele_2_col])
            normal_allele_3 = safe_convert_to_int(split_line[allele_3_col])
            if normal_allele_3 != 0:
                continue
            else:
                if normal_allele_2 == 0:
                    normal_allele_2 = normal_allele_1 # if homozygous, count the allele twice
                for allele in [normal_allele_1, normal_allele_2]:
                    if allele != 0:
                        matrix[line_num, allele]+=1
    e = time.time()
    print(e-st)

    if (line_num+1)!=LOCI_LINE_COUNT:
        raise RuntimeError(f"Only had {line_num} lines. Should have {LOCI_LINE_COUNT}")
    return matrix

def create_loci_matrix():
    matrix = np.zeros((LOCI_LINE_COUNT, 41), dtype=np.int32)
    return matrix


def create_matrices_for_tsvs(msmutect_output_files: List[str]) -> np.ndarray:
    loci_matrix = create_loci_matrix()
    for f in msmutect_output_files:
        loci_matrix = create_matrix_for_tsv(f, loci_matrix)
    return loci_matrix


def list_of_gz_files(src_dir: str):
    return glob.glob(os.path.join(src_dir, "*.full.mut.tsv.gz"))


def main(src_dir: str, start_idx: int, stop_idx: int):
    msmutect_files = list_of_gz_files(src_dir)[start_idx:stop_idx]
    print(msmutect_files)
    ret = create_matrices_for_tsvs(msmutect_files)
    np.save(f"{os.path.basename(src_dir)}_{start_idx}_{stop_idx}.npy", ret)


def main2(input_file: str, output_file: str):
    ret = create_matrices_for_tsvs([input_file])
    np.save(output_file, ret)


if __name__ == '__main__':
    main2(sys.argv[1], sys.argv[2])
    #main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]))