import csv, time
from typing import List, Union, Dict

from scipy.stats import binom

from positives_simulation.COLUMN_NAMES import COLUMN_NAMES
from positives_simulation.DBwriteRequest import DBwriteRequest
from positives_simulation.Distribution import Distribution
from positives_simulation.RandomAllelesDB import RandomAllelesDB


def is_round_number(number: str) -> bool:
    number_float = float(number)
    number_int = int(number_float)
    return abs(number_float-number_int) < 0.001 # for rounding point error

def return_default_if_na_else_return(x: str, default) -> Union[int, None]:
    if x == "NA":
        return default
    else:
        return int(x)


def cdf_test(first_allele_reads: int, second_allele_reads: int, p_equal: float = 0.3):
    return  binom.cdf(min(first_allele_reads, second_allele_reads), first_allele_reads + second_allele_reads, 0.5) > p_equal


def reads_equally_distributed(reads_1: str, reads_2: str) -> bool:
    if "NA" in (reads_1, reads_2):
        return True
    else:
        return cdf_test(int(reads_1), int(reads_2))



def repeat_lengths_from_row(row: dict) -> Dict[int, int]:
    column_names = [("NORMAL_MOTIF_REPEATS_1", "NORMAL_SUPPORTING_READS_1"),
                    ("NORMAL_MOTIF_REPEATS_2", "NORMAL_SUPPORTING_READS_2"),
                    ("NORMAL_MOTIF_REPEATS_3", "NORMAL_SUPPORTING_READS_3"),
                    ("NORMAL_MOTIF_REPEATS_4", "NORMAL_SUPPORTING_READS_4"),
                    ("NORMAL_MOTIF_REPEATS_5", "NORMAL_SUPPORTING_READS_5"),
                    ("NORMAL_MOTIF_REPEATS_6", "NORMAL_SUPPORTING_READS_6")]
    repeat_lengths = dict()
    for pair in column_names:
        if row[pair[0]] == "NA":
            return repeat_lengths
        else:
            repeat_lengths[int(row[pair[0]])] = int(row[pair[1]])
    return repeat_lengths

def flexint(x: str) -> int:
    return int(float(x))

def create_distribution_from_row(row: dict) -> Distribution:
    return Distribution(row["PATTERN"], int(row["NORMAL_ALLELE_1"]), return_default_if_na_else_return(row["NORMAL_ALLELE_2"], None),
                        int(row["END"])-int(row["START"])+1, flexint(row["REFERENCE_REPEATS"]), repeat_lengths_from_row(row))


def create_write_request_from_row(row: dict) -> Union[DBwriteRequest, None]:
    if ((row["NORMAL_ALLELE_1"] == "NA") or (row["NORMAL_ALLELE_3"] != "NA") or
            (not is_round_number(row["REFERENCE_REPEATS"])) or (len(row["PATTERN"])>1)) or \
                (not reads_equally_distributed(row["NORMAL_SUPPORTING_READS_1"], row["NORMAL_SUPPORTING_READS_2"])): # must have 1-2 alleles
        return None
    heterozygous = (row["NORMAL_ALLELE_2"] != "NA")
    is_reference = flexint(row["REFERENCE_REPEATS"]) in [int(row["NORMAL_ALLELE_1"]),
                                                         int(return_default_if_na_else_return(row["NORMAL_ALLELE_2"],
                                                                                              0))]
    if heterozygous:
        keys = [int(row["NORMAL_ALLELE_1"]), int(row["NORMAL_ALLELE_2"])]
    else:
        keys = [int(row["NORMAL_ALLELE_1"])]
    distribution = create_distribution_from_row(row)
    return DBwriteRequest(keys, not heterozygous, is_reference, distribution)


def create_patient_write_requests(patient_filepath: str) -> List[DBwriteRequest]:
    write_requests = []
    with open(patient_filepath, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        st = time.time()
        for i, row in enumerate(reader):
            # print(row["REFERENCE_REPEATS"])
            # print(row)
            new_write_request = create_write_request_from_row(row)

            if new_write_request is not None:
                write_requests.append(new_write_request)
            if i % 100_000 == 0:
                e = time.time()
                print(e-st)
                st=time.time()
            if i == 500_000:
                return write_requests
    return write_requests



def add_patient_to_db(patient_name: str, input_patient_path: str):
    alleles_db = RandomAllelesDB("alleles.db")
    write_requests = create_patient_write_requests(input_patient_path)
    alleles_db.add_patient(patient_name, write_requests)
    x=1

def read_db(patient_name: str):
    alleles_db = RandomAllelesDB("alleles.db")

if __name__ == "__main__":
    read_db("hg001")
    # add_patient_to_db("hg001", "/home/avraham/MaruvkaLab/msmutect_postprocessing/data/full_gib_files/msmutect_normal0_tumor0.full.mut.tsv")
