import random
from typing import Dict, Union, List
from dataclasses import dataclass, field

from positives_simulation.COLUMN_NAMES import COLUMN_NAMES
from positives_simulation.repeat_support_handling import dict_to_csv_representation, return_na_if_none_else_return_input


@dataclass
class Distribution:
    pattern: str
    allele_1: int
    allele_2: Union[int, None]
    locus_length: int
    reference_size: int
    repeat_lengths: Dict[int, int]
    _probability_map: List[int] = field(default_factory=list)

    def num_reads(self):
        return sum(self.repeat_lengths.values())

    def construct_probability_map(self):
        probability_map = []
        for motif_repeats in self.repeat_lengths:
            for weight in range(self.repeat_lengths[motif_repeats]):
                probability_map.append(motif_repeats)
        self._probability_map = probability_map

    def randomly_select_read(self):
        if self._probability_map is None or ((self._probability_map) != sum(self.repeat_lengths.values())):
            self.construct_probability_map()
        return self._probability_map[random.randint(0, len(self._probability_map) - 1)]

    def __str__(self):
        return f"{self.pattern},{self.allele_1},{return_na_if_none_else_return_input(self.allele_2)},{self.locus_length},{self.reference_size},{dict_to_csv_representation(self.repeat_lengths)}"

    @staticmethod
    def column_names():
        motif_repeats_columns = [COLUMN_NAMES.REPEAT_LENGTH_1,
                                 COLUMN_NAMES.REPEAT_LENGTH_2,
                                 COLUMN_NAMES.REPEAT_LENGTH_3,
                                 COLUMN_NAMES.REPEAT_LENGTH_4,
                                 COLUMN_NAMES.REPEAT_LENGTH_5,
                                 COLUMN_NAMES.REPEAT_LENGTH_6]
        motif_repeat_support_columns = [
            COLUMN_NAMES.REPEAT_SUPPORT_1,
            COLUMN_NAMES.REPEAT_SUPPORT_2,
            COLUMN_NAMES.REPEAT_SUPPORT_3,
            COLUMN_NAMES.REPEAT_SUPPORT_4,
            COLUMN_NAMES.REPEAT_SUPPORT_5,
            COLUMN_NAMES.REPEAT_SUPPORT_6,
        ]

        return "PATTERN,ALLELE_1,ALLELE_2,LOCUS_LENGTH,REFERENCE_NUM_REPEATS,"+",".join(motif_repeats_columns)+","+",".join(motif_repeat_support_columns)

if __name__ == "__main__":
    print(Distribution.column_names())