import os, re, csv, random, sys
from collections import defaultdict
from dataclasses import dataclass
from typing import List, Dict

from positives_simulation.COLUMN_NAMES import COLUMN_NAMES
from positives_simulation.DBwriteRequest import DBwriteRequest
from positives_simulation.Distribution import Distribution


@dataclass
class PatientDBengine:
    def __init__(self, src_dir: str, create=False) -> None:
        self.db = dict()
        self.src_dir = src_dir
        all_files = os.listdir(self.src_dir)
        # should make create mode write only, visa versa for create=False
        if not create:
            if len(all_files) == 0 :
                raise RuntimeError(f"No files found")
            self.num_keys = all_files[0].count("_")+1
            self.validate_files(all_files, num_keys=self.num_keys)
            for f in all_files:
                self.db[f[:f.find(".csv")]] = self.load_file(os.path.join(src_dir, f))

    def has_entry(self, keys: List[int]) -> bool:
        self.validate_num_keys(keys)
        relevant_distributions = self.db["_".join([str(key) for key in keys])]
        return len(relevant_distributions) > 0

    def get_random_entry(self, keys: List[int]) -> Distribution:
        self.validate_num_keys(keys)
        relevant_distributions = self.db["_".join([str(key) for key in keys])]
        if len(relevant_distributions) == 0:
            raise RuntimeError(f"No distributions found for given keys")
        random_num = random.randint(0,len(relevant_distributions)-1)
        return relevant_distributions[random_num]

    def validate_num_keys(self, keys: List[int]) -> None:
        if len(keys) != self.num_keys:
            raise RuntimeError(f"Number of given keys does not match number of keys")


    def construct_repeat_lengths_dict(self, row) -> Dict[int, int]:
        column_names = [(COLUMN_NAMES.REPEAT_LENGTH_1, COLUMN_NAMES.REPEAT_SUPPORT_1),
                        (COLUMN_NAMES.REPEAT_LENGTH_2, COLUMN_NAMES.REPEAT_SUPPORT_2),
                        (COLUMN_NAMES.REPEAT_LENGTH_3, COLUMN_NAMES.REPEAT_SUPPORT_3),
                        (COLUMN_NAMES.REPEAT_LENGTH_4, COLUMN_NAMES.REPEAT_SUPPORT_4),
                        (COLUMN_NAMES.REPEAT_LENGTH_5, COLUMN_NAMES.REPEAT_SUPPORT_5),
                        (COLUMN_NAMES.REPEAT_LENGTH_6, COLUMN_NAMES.REPEAT_SUPPORT_6)]
        repeat_lengths = dict()
        for pair in column_names:
            if row[pair[0]] == "NA":
                return repeat_lengths
            else:
                repeat_lengths[int(pair[0])] = int(row[pair[1]])
        return repeat_lengths

    def zero_if_NA_otherwise_self(self, x: str) -> str:
        return "0" if x == "NA" else x

    def load_file(self, filepath: str) -> List[Distribution]:
        distributions = []
        with open(filepath, "r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                distributions.append(Distribution(row[COLUMN_NAMES.PATTERN],
                                                  int(row[COLUMN_NAMES.ALLELE_1]),
                                                  int(self.zero_if_NA_otherwise_self(row[COLUMN_NAMES.ALLELE_2])),
                                                  int(row[COLUMN_NAMES.LOCUS_LENGTH]),
                                                  int(row[COLUMN_NAMES.REFERENCE]),
                                                  self.construct_repeat_lengths_dict(row)))
        return distributions

    def validate_files(self, file_list: List[str], num_keys: int) -> None:
        pattern = "^" + "_".join(["\d+" for _ in range(num_keys)]) + "\.csv$"
        compiled_pattern = re.compile(pattern)
        for f in file_list:
            if not compiled_pattern.match(f):
                raise RuntimeError(f"File {f} does not match pattern {pattern}. Very likely has incorrect number of keys")

    def size(self):
        return sum([sys.getsizeof(lst) for lst in self.db.values()])

    def write_entries(self, write_requests: List[DBwriteRequest]) -> None:
        db = defaultdict(list)
        for write_request in write_requests:
            db["_".join([str(x) for x in write_request.keys])].append(str(write_request.distribution))

        for key_set in db:
            with open(os.path.join(self.src_dir, f"{key_set}.csv"), "w+") as opened_file:
                opened_file.write(Distribution.column_names()+"\n")
                opened_file.write("\n".join(db[key_set]))






if __name__ == '__main__':
    a = PatientDBengine("tmptest")
