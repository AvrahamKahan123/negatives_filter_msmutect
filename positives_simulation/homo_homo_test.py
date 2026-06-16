#!/Local/bfe_maruvka/anaconda3/bin/python -u
import random, time

import numpy as np
from typing import List, Dict
import sys

from collections import defaultdict

from positives_simulation.MSMuTect_4.src.IndelCalling.Locus import Locus

# sys.path.append("/home/avraham/MaruvkaLab/msmutect_postprocessing/positives_simulation/MSMuTect_4") # hideous, short term solution
sys.path.append("/storage/bfe_maruvka/avrahamk/Positives/msmutect_postprocessing/positives_simulation/MSMuTect_4") # hideous, short term solution

from positives_simulation.MSMuTect_4.src.GenomicUtils.NoiseTable import get_noise_table
from positives_simulation.MSMuTect_4.src.IndelCalling import CallAllelesFast, CallMutations
from positives_simulation.MSMuTect_4.src.IndelCalling.Histogram import Histogram
from positives_simulation.MSMuTect_4.src.IndelCalling.MutationCall import MutationCall

from positives_simulation.Distribution import Distribution
from positives_simulation.RandomAllelesDB import RandomAllelesDB
from positives_simulation.RandomRequest import RandomRequest


def all_purities() -> List[float]:
    purity_increment = 0.05
    return [round(x, 2) for x in np.arange(purity_increment, 1 + purity_increment, purity_increment)]

ALL_PURITIES = all_purities()
NOISE_TABLE = get_noise_table()
MINIMUM_NUM_REPEATS=5
MAXIMUM_NUM_REPEATS=40
REQUIRED_READ_SUPPORT=5

def initiate_matrix() -> np.ndarray:
    return np.zeros((MAXIMUM_NUM_REPEATS, MAXIMUM_NUM_REPEATS), dtype=np.int32)

def move_probabilities_to_absolute(probabilities: List[float]):
    current = 0
    absolute_probabilities = []
    for p in probabilities:
        absolute_probabilities.append(p+current)
        current += p
    return absolute_probabilities

def simulate_distribution(distributions: List[Distribution], fractions: List[float], num_reads: int) -> Dict[int, int]:
    if abs(sum(fractions) - 1) > 1e-3:
        raise RuntimeError("Probabilities do not sum to 1")
    adjusted_probabilities = move_probabilities_to_absolute(fractions)
    simulated_distribution = defaultdict(int)
    for _ in range(num_reads):
        p = random.random()
        found_dist = False
        for adjusted_probability, distribution in zip(adjusted_probabilities, distributions):
            if p < adjusted_probability:
                randomly_selected_length = distribution.randomly_select_read()
                simulated_distribution[randomly_selected_length] += 1
                found_dist = True
                break
        if not found_dist:
            raise RuntimeError("No distribution found. Probabilities do not sum to 1")
    return dict(simulated_distribution)


def construct_histogram_from_distribution(distribution: Dict[int, int], reference_size: int) -> Histogram:
    locus = Locus("1", 1, reference_size, "A"*reference_size, reference_size, sequence="A"*reference_size) # fake locus
    histogram = Histogram(locus)
    histogram.repeat_lengths = defaultdict(int)
    for length, count in distribution.items():
        histogram.repeat_lengths[length] = count
    return histogram


def is_mutation(normal_distribution: Dict[int, int], tumor_distribution: Dict[int, int], reference_size: int) -> bool:
    normal_histogram = construct_histogram_from_distribution(normal_distribution, reference_size)
    tumor_histogram = construct_histogram_from_distribution(tumor_distribution, reference_size)
    tumor_alleles = CallAllelesFast.calculate_alleles(tumor_histogram, NOISE_TABLE,
                                                      required_read_support=REQUIRED_READ_SUPPORT)
    normal_alleles = CallAllelesFast.calculate_alleles(normal_histogram, NOISE_TABLE,
                                                       required_read_support=REQUIRED_READ_SUPPORT)
    mutation_call = CallMutations.call_mutations(normal_alleles, tumor_alleles, NOISE_TABLE)
    return mutation_call.call == MutationCall.MUTATION


def find_neccesary_read_support(reference_distributions: List[Distribution],
                                non_reference_distributions: List[Distribution],
                                purity: float, tumor_factor: float = 2.5, mutation_percent_threshold: float = 0.5) -> int:
    if len(non_reference_distributions) != len(reference_distributions):
        raise ValueError('The number of non_reference_distributions must match the number of reference_distributions')
    tumor_fraction = purity / 2
    normal_fraction = 1-tumor_fraction
    num_simulations = len(reference_distributions)
    for num_reads in range(5, 150):
        num_mutations = 0
        for reference_distribution, non_reference_distribution in zip(reference_distributions, non_reference_distributions):
            randomly_simulated_normal = simulate_distribution([reference_distribution], [1.0],
                                                            num_reads)
            randomly_simulated_tumor = simulate_distribution([reference_distribution, non_reference_distribution],
                                                      [normal_fraction, tumor_fraction], int(num_reads*tumor_factor))
            if is_mutation(randomly_simulated_normal, randomly_simulated_tumor, reference_size=reference_distribution.reference_size):
                num_mutations+=1
        if ((num_mutations/num_simulations)+0.0001) >= mutation_percent_threshold: # so we won't have floating point problem
            return num_reads
    return 1e4

def get_n_random_distributions(n: int, alleles_db: RandomAllelesDB, patient_name: str, homozygous: bool,
                               reference_size: int, non_reference_size: int = None) -> List[Distribution]:
    if non_reference_size is None:
        request = RandomRequest(patient_name, homozygous, reference_size)
    else:
        request = RandomRequest(patient_name, homozygous, reference_size, non_reference_size)
    test_response = alleles_db.get(request)
    if not test_response.succeeded:
        return []
    else:
        distributions = [answer.distribution for answer in [alleles_db.get(request) for _ in range(n)]]
        return distributions

def simulate_purity(alleles_db: RandomAllelesDB, patient_name: str, homozygous: bool,
                    purity: float, num_simulations: int) -> np.ndarray:

    matrix = initiate_matrix()
    for reference_size in range(MINIMUM_NUM_REPEATS, MAXIMUM_NUM_REPEATS+1):
        for mutation_size in range(MINIMUM_NUM_REPEATS, MAXIMUM_NUM_REPEATS+1):
            reference_distributions = get_n_random_distributions(num_simulations, alleles_db, patient_name,
                                                                 homozygous, reference_size)
            non_reference_distributions = get_n_random_distributions(num_simulations, alleles_db, patient_name,
                                                                 homozygous, reference_size, mutation_size)
            if len(non_reference_distributions) != 0 and len(reference_distributions) != 0:
                st = time.time()
                matrix[reference_size-1, mutation_size-1] = (
                    find_neccesary_read_support(reference_distributions, non_reference_distributions, purity))
                e=time.time()
                print(f"Simulation took {e-st} seconds")
    return matrix

def main():
    num_simulations = 1
    patient_name = "hg001"
    homozygous = True
    alleles_db = RandomAllelesDB("alleles.db")
    alleles_db.load_patient(patient_name)
    for purity in ALL_PURITIES:
        current_purity_matrix = simulate_purity(alleles_db, patient_name, homozygous, purity, num_simulations)
        np.save(f"purity={purity}.npy", current_purity_matrix)
        print(f"Saved purity {purity}")

def main_parralel():
    num_simulations = 50
    patient_name = "hg001"
    homozygous = True
    alleles_db = RandomAllelesDB("alleles.db")
    alleles_db.load_patient(patient_name)
    purity = float(sys.argv[1])
    print(purity)
    current_purity_matrix = simulate_purity(alleles_db, patient_name, homozygous, purity, num_simulations)
    np.save(f"purity={purity}.npy", current_purity_matrix)
    print(f"Saved purity {purity}")

if __name__ == '__main__':
    main_parralel()