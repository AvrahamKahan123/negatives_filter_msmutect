import numpy as np
from typing import List, Dict

from positives_simulation.Distribution import Distribution
from positives_simulation.RandomAllelesDB import RandomAllelesDB


def all_purities() -> List[float]:
    purity_increment = 0.05
    return [round(x, 2) for x in np.arange(purity_increment, 1 + purity_increment, purity_increment)]

def initiate_matrix() -> np.ndarray:
    return np.zeros((40, 40), dtype=np.int32)

def simulate_distribution(distributions: List[Distribution], fractions: List[float], num_reads: int) -> Dict[int, int]:
    pass

def is_mutation() -> bool:
    raise NotImplementedError


def find_neccesary_read_support(reference_distribution: Distribution, non_reference_distribution: Distribution,
                                purity: float, tumor_factor: float = 2.5, num_simulations = 20,
                                mutation_percent_threshold: float = 0.5) -> int:
    tumor_fraction = purity / 2
    normal_fraction = 1-tumor_fraction
    for num_reads in range(5, 150):
        num_mutations = 0
        for simulation in range(num_simulations):
            randomly_simulated_normal = simulate_distribution([reference_distribution], [normal_fraction],
                                                            num_reads)
            randomly_simulated_tumor = simulate_distribution([reference_distribution, non_reference_distribution],
                                                      [normal_fraction, tumor_fraction], int(num_reads*tumor_factor))
            if is_mutation(randomly_simulated_normal, randomly_simulated_tumor):
                num_mutations+=1
        if (num_mutations/num_simulations) > mutation_percent_threshold:
            return num_mutations
    return np.inf


def main():
    alleles_db = RandomAllelesDB("alleles.db")
    alleles_db.load_patient("hg001")


if __name__ == '__main__':
    main()