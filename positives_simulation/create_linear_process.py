#!/Local/bfe_maruvka/anaconda3/bin/python -u
import random, time
from multiprocessing import Pool

import numpy as np
from typing import List, Dict
import sys

from collections import defaultdict

from positives_simulation.MSMuTect_4.src.IndelCalling.Locus import Locus

sys.path.append("/home/avraham/MaruvkaLab/msmutect_postprocessing/positives_simulation/MSMuTect_4") # hideous, short term solution
# sys.path.append("/storage/bfe_maruvka/avrahamk/Positives/msmutect_postprocessing/positives_simulation/MSMuTect_4") # hideous, short term solution

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

REQUIRED_READ_SUPPORT=5
ALL_PURITIES = all_purities()
NOISE_TABLE = get_noise_table()

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





def get_random_reference_length():
    probabilities = [(5, 0.6275589627651392), (6, 0.8057150885640911), (7, 0.8727316026797027), (8, 0.8958289689765649),
                     (9, 0.9085530069746692), (10, 0.9411773215570167), (11, 0.9512084903891368), (12, 0.9596388253666424),
                     (13, 0.9667088086514611), (14, 0.9717296806536059), (15, 0.9757163374415663), (16, 0.97970542056773),
                     (17, 0.9825336610934745), (18, 0.9850729842531565), (19, 0.9873540168440975), (20, 0.9892791101986044),
                     (21, 0.9909351089123705), (22, 0.9924950154772908), (23, 0.9937975012928261), (24, 0.994932295092534),
                     (25, 0.9959213254801768), (26, 0.99670265216162), (27, 0.9973368420320151), (28, 0.9978811200489881),
                     (29, 0.9983066448338955), (30, 0.9986437685042622), (31, 0.9989072779891387), (32, 0.9991164466542535),
                     (33, 0.9992897696059113), (34, 0.9994411181738407), (35, 0.9995664179410587), (36, 0.9996743671011246),
                     (37, 0.9997777840445464), (38, 0.9998655442395618), (39, 0.9999394788848146), (40, 0.9999999999999999)]

    rand_float = random.random()
    for p in reversed(probabilities):
        if rand_float > p[1]:
            return p[0]
    return probabilities[0][0]

def get_random_mutation_length(p_mut_map: np.ndarray, reference_length: int) -> int:
    relevant_row = p_mut_map[reference_length]
    rand_float = random.random()
    for i, p in enumerate(reversed(relevant_row)):
        if rand_float > p:
            return len(relevant_row) - i
    raise RuntimeError("Should never get here")

def purity_mutation_rate(purity: float, num_simulation: int = 500) -> float:
    p_mut_map = np.load("/home/avraham/MaruvkaLab/msmutect_postprocessing/positives_simulation/p_mut_map.npy")
    tumor_fraction = purity/2
    normal_fraction = 1-tumor_fraction
    tumor_factor = 2.5
    alleles_db = RandomAllelesDB("alleles.db")
    alleles_db.load_patient("hg001")
    num_mutations = 0
    # for sim in range(num_simulation):
    sim = 0
    while sim < num_simulation:
        reference_length = get_random_reference_length()
        random_mutation_length = get_random_mutation_length(p_mut_map, reference_length)
        normal_request = RandomRequest("hg001", True, reference_length)
        normal_response = alleles_db.get(normal_request)
        normal_distribution = normal_response.distribution
        tumor_request = RandomRequest("hg001", True, reference_length, random_mutation_length)
        tumor_response = alleles_db.get(tumor_request)
        if not tumor_response.succeeded or not normal_response.succeeded:
            # print(f"FAILED DB QUERY: REF: {reference_length}, MUT: {random_mutation_length}")
            continue
        else:
            sim+=1
        tumor_distribution = tumor_response.distribution
        num_reads = normal_distribution.num_reads()
        randomly_simulated_normal = simulate_distribution([normal_distribution], [1.0],
                                                          num_reads)
        randomly_simulated_tumor = simulate_distribution([normal_distribution, tumor_distribution],
                                                         [normal_fraction, tumor_fraction],
                                                         int(num_reads * tumor_factor))
        if is_mutation(randomly_simulated_normal, randomly_simulated_tumor,
                       reference_size=reference_length):
            num_mutations += 1

    return num_mutations/num_simulation


def main():

    m_rates = []
    for purity in ALL_PURITIES:
        mutation_rate = purity_mutation_rate(purity)
        m_rates.append(mutation_rate)
        # print(mutation_rate)
    print(m_rates)

def call_main(_):
    # ignore the argument, just call myfunc
    return main()

if __name__ == '__main__':
    st = time.time()
    with Pool(16) as p:
        results = p.map(call_main, range(16))
    e = time.time()
    print(e-st)