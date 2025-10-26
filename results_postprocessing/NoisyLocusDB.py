import os
import numpy as np
from results_postprocessing.cancer_data import data_directory


def db_filepath():
    return os.path.join(data_directory(), "noisy_loci.txt")


class NoisyLocusDB:
    def __init__(self):
        self.db = set()
        with open(db_filepath(), 'r') as db_file:
            for line in db_file:
                self.db.add(line.rstrip())

    def locus_present_in_db(self, chromosome: str, start: int, stop: int, pattern: str):
        encoded_query = self.encode(chromosome, start, stop, pattern)
        return encoded_query in self.db

    @staticmethod
    def encode(chromosome: str, start: int, stop: int, pattern: str):
        return f"{chromosome}_{start}_{stop}_{pattern}"


def create_db(loci_fp: str, noisy_locus_npy: str):
    noisy_loci_ref = np.load(noisy_locus_npy)
    db = open(db_filepath(), 'w+')
    with (open(loci_fp, 'r') as loci_file):
        for line_num, line in enumerate(loci_file):
            if noisy_loci_ref[line_num]:
                split_line = line.split("\t")
                db.write(NoisyLocusDB.encode(split_line[0].replace("chr", ""), split_line[3], split_line[4], split_line[12]))
                db.write("\n")
    db.close()


if __name__ == '__main__':
    create_db("C:\\Users\\avrah\MaruvkaLab\msmutect_runs\data\GRCh38.d1.vd1_1to15_repetitive_loci_sorted_fixed.phobos",
              os.path.join(data_directory(), "noisy_locus.npy"))
    # v=NoisyLocusDB()
