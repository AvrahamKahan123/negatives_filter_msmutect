import random, time, sys
from dataclasses import dataclass
from typing import List, Dict, Tuple
import numpy as np


class HistogramSet:
    def __init__(self, alleles: List[int] = None, repeat_lengths: Dict[int, int] = None):
        # alleles are the alleles called by the algorithm
        # repeat_lengths is of form {MOTIF: MOTIF_SUPPORT}
        self.alleles = alleles
        if repeat_lengths is None:
            self._repeat_lengths = dict()
        else:
            self._repeat_lengths = repeat_lengths
        self._probability_map: np.ndarray = None # is array of form [REPEAT_LENGTH_1 (REPEAT_LENGTH_1_SUPPORT times), REPEAT_LENGTH_2 (REPEAT_LENGTH_2_SUPPORT times), ...]
        self._sorted_repeat_lengths: Dict[int, int] = None

    def homozygous(self):
        if type(self.alleles) is None:
            raise RuntimeError("No alleles set")
        return len(self.alleles) == 1

    def remove_repeat_length(self, length: int):
        if length in self._repeat_lengths:
            del self._repeat_lengths[length]

    def num_reads(self) -> int:
        return sum(list(self._repeat_lengths.values()))

    def add_read(self, repeat_length: int):
        if repeat_length in self._repeat_lengths:
            self._repeat_lengths[repeat_length]+=1
        else:
            self._repeat_lengths[repeat_length]=1
        self._probability_map = None # is now invalid

    def create_probability_map(self) -> None:
        pmap = np.zeros((self.num_reads()), dtype=np.int32)
        current_index = 0
        for repeat, repeat_support in self._repeat_lengths.items():
            pmap[current_index:current_index+repeat_support] = repeat
            current_index+=repeat_support
        self._probability_map = pmap

    def randomly_select_read(self, probability: float = None) -> int:
        if self._probability_map is None:
            self.create_probability_map()
        if self.num_reads()==0:
            raise RuntimeError("Cannot randomly select from an empty histogram")
        random_int = random.randint(0, self.num_reads()-1)
        return self._probability_map[random_int]

    def copy(self):
        return HistogramSet(alleles=self.alleles, repeat_lengths=self._repeat_lengths.copy())

    def create_sorted_repeat_lengths(self):
        sorted_keys = list(sorted(self._repeat_lengths, key=self._repeat_lengths.get, reverse=True))
        ret = {k: self._repeat_lengths[k] for k in sorted_keys}
        return ret


    @property
    def repeat_lengths(self) -> Dict[int, int]:
        if self._sorted_repeat_lengths is None:
            self._sorted_repeat_lengths = self.create_sorted_repeat_lengths()
        return self._sorted_repeat_lengths



@dataclass
class Result:
    result: object
    succeeded: bool
    reason: str


def probabilities_to_map(probabilities: List[float]) -> List[float]:
    absolute_probabilities = [probabilities[0]]
    for p in probabilities[1:]:
        absolute_probabilities.append(p+absolute_probabilities[-1])
    return absolute_probabilities


def mix_histograms_probabilistically(histogram_sets: List[HistogramSet], relative_probabilities: List[float], read_support: int) -> HistogramSet:
    if abs(sum(relative_probabilities)-1) > 1e-5:
        raise RuntimeError("Probabilities do not sum to 1")
    probability_map = probabilities_to_map(relative_probabilities)
    mixed_set = HistogramSet()
    for _ in range(read_support):
        p = random.random()
        new_read = None
        for hist_set, prob in zip(histogram_sets, probability_map):
            if p < prob:
                new_read = hist_set.randomly_select_read()
                break
        if new_read is None:
            raise RuntimeError(f"Probability {p} did not fit map {probability_map}. Most likely cause is a rounding error")
        else:
            mixed_set.add_read(new_read)
    return mixed_set


def find_less_dominant_normal_allele(normal: HistogramSet) -> int:
    return list(normal.repeat_lengths.keys())[1] # is normal, so can only have 2 alleles


def mix_histograms(normal: HistogramSet, tumor: HistogramSet, purity: float) -> HistogramSet:
    # will return None if the normal is heterozygous but could not figure out predominant repeat lengh
    read_support = tumor.num_reads()
    true_purity = purity/2 # we divide purity by 2 because we assume the mutation is heterozygous
    if normal.homozygous():
        normal_allele = normal.alleles[0]
        tumor_with_normal_allele_removed = tumor.copy()
        tumor_with_normal_allele_removed.remove_repeat_length(normal_allele)
        hist = mix_histograms_probabilistically([normal, tumor_with_normal_allele_removed], [1-true_purity, true_purity], read_support)
    else:
        allele_replaced_with_mutation = find_less_dominant_normal_allele(normal)
        normal_without_less_dominant = normal.copy()
        normal_without_less_dominant.remove_repeat_length(allele_replaced_with_mutation)
        tumor_with_normal_allele_removed = tumor.copy()
        tumor_with_normal_allele_removed.remove_repeat_length(allele_replaced_with_mutation)
        hist = mix_histograms_probabilistically([normal, normal_without_less_dominant,
                                                 tumor_with_normal_allele_removed], [1-purity, true_purity, true_purity], read_support)
    return hist


def convert_str_list_to_ints(lst: List[str]):
    return [int(x) for x in lst if x != "0"]


def histogram_sets(line: str) -> Tuple[HistogramSet, HistogramSet]:
    split_line = line.split("\t")
    normal_alleles = convert_str_list_to_ints(split_line[20:24])
    normal_repeats = convert_str_list_to_ints(split_line[6:12])
    normal_repeat_support = convert_str_list_to_ints(split_line[12:18])
    tumor_repeats = convert_str_list_to_ints(split_line[28:34])
    tumor_repeats_support = convert_str_list_to_ints(split_line[34:40])
    normal_histogram_set = HistogramSet(normal_alleles, {repeat: repeat_support for repeat, repeat_support
                                                         in zip(normal_repeats, normal_repeat_support)})
    tumor_histogram_set = HistogramSet(repeat_lengths={repeat: repeat_support for repeat, repeat_support
                                                         in zip(tumor_repeats, tumor_repeats_support)})
    return normal_histogram_set, tumor_histogram_set


class MemoryEfficientFileWriter:
    def __init__(self, output_fp: str, line_threshold: int = 100_000):
        self.output_file = open(output_fp, 'w+')
        self.line_threshold = line_threshold
        self.lines_queue = []

    def writeline(self, line: str):
        self.lines_queue.append(line)
        if len(self.lines_queue) > self.line_threshold:
            self.flush_queue()

    def flush_queue(self):
        self.output_file.write("\n".join(self.lines_queue)+"\n")
        self.lines_queue = []

    def __del__(self):
        self.flush_queue()
        self.output_file.close()


def get_first_n_entries(line: str, entry_idx: int):
    return "\t".join(line.split("\t")[:entry_idx])


def format_header_line(header_line: str):
    return get_first_n_entries(header_line, 18)


def fill_in_lists(lst: list, desired_length: int, fill_in_val: object) -> list:
    missing_length = max(desired_length - len(lst), 0) # at least 0
    ret = lst + [fill_in_val for _ in range(missing_length)]
    return ret[:desired_length]


def format_histogram_set(histogram_set: HistogramSet):
    lengths = fill_in_lists(list(histogram_set.repeat_lengths.keys()), 6, "NA")
    support = fill_in_lists(list(histogram_set.repeat_lengths.values()), 6, "NA")
    return "\t".join([str(x) for x in lengths+support]) # append the lists


def create_mixed_histograms_file(input_fp: str, output_fp: str, purity: float, number_of_entries: int):
    output_file = MemoryEfficientFileWriter(output_fp)
    with open(input_fp, 'r') as opened_mut_file:
        header_line = opened_mut_file.readline()
        output_file.writeline(format_header_line(header_line))
        for line in opened_mut_file:
            normal_set, tumor_set = histogram_sets(line)
            line_start = get_first_n_entries(line, 6) +"\t"
            for _ in range(number_of_entries):
                new_tumor_set = mix_histograms(normal_set, tumor_set, purity=purity)
                formatted_tumor_set = format_histogram_set(new_tumor_set)
                output_file.writeline(line_start+formatted_tumor_set)


def create_randomized_purity_files(input_fp: str, output_prefix: str, number_of_entries_per_purity: int):
    step = 0.05
    for purity in [round(x, 2) for x in np.arange(0.1, 1+step, step)]:
        create_mixed_histograms_file(input_fp,
                                     f"{output_prefix}_{number_of_entries_per_purity}x_{purity}purity.hist.tsv",
                                     purity, number_of_entries_per_purity)


if __name__ == '__main__':
    # main("../data/msi_example_files/TCGA-A6-5661.called.filt.mut.tsv",
    #      "TCGA-A6-5661", 10)
    st = time.time()
    create_randomized_purity_files("/home/avraham/MaruvkaLab/msmutect_postprocessing/data/msi_example_files/TCGA-A6-5661.called.filt.mut.tsv",
                                   "/home/avraham/MaruvkaLab/msmutect_postprocessing/data/positives/TCGA-A6-5661.rand", 1)

    # create_randomized_purity_files("/home/avraham/MaruvkaLab/msmutect_postprocessing/data/mss_example_files/TCGA-2G-AALO.called.filt.mut.tsv",
    #                                "/home/avraham/MaruvkaLab/msmutect_postprocessing/data/positives/TCGA-2G-AALO.rand", 1)

    # create_randomized_purity_files(sys.argv[1], sys.argv[2], 1)
    e = time.time()
    print(e-st)

