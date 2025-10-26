import os
from dataclasses import dataclass, field
from typing import List, Dict

from results_postprocessing.enums import MSI_CLASSIFICATION


@dataclass
class Sample:
    mutations_tsv_fp: str
    cancer_type: str
    classification: str

    def patient_name(self):
        # differs from sample_name in that it just returns the name of the patient; multiple
        return self.patient_name_from_filename(self.mutations_tsv_fp)

    def sample_name(self):
        return self.get_filename_no_ext(self.mutations_tsv_fp)

    @staticmethod
    def patient_name_from_filename(name: str):
        suffix_idx = Sample.find_first_period_or_underscore(name)
        if suffix_idx == -1: # no period/underscore
            return name
        else:
            return name[:suffix_idx]

    @staticmethod
    def find_first_period_or_underscore(s: str) -> int:
        """
        Find the index of the first period (.) or underscore (_) in a string.
        """
        # Find index of '.' and '_'
        period_index = s.find('.')
        underscore_index = s.find('_')

        if period_index == -1:
            return underscore_index
        if underscore_index == -1:
            return period_index
        return min(period_index, underscore_index)

    @staticmethod
    def get_filename_no_ext(fp: str):
        fp = os.path.basename(fp)
        period_loc = fp.find(".")
        return fp[:period_loc]


@dataclass
class SamplesSet:
    cancer_type: str
    mss_patients: List[Sample] = field(default_factory=list)
    msi_patients: List[Sample] = field(default_factory=list)
    negative_patients: List[Sample] = field(default_factory=list)

    def add_patient(self, sample: Sample):
        if sample.classification == MSI_CLASSIFICATION.MSI:
            self.msi_patients.append(sample)
        elif sample.classification == MSI_CLASSIFICATION.MSS:
            self.mss_patients.append(sample)
        elif sample.classification == MSI_CLASSIFICATION.NEGATIVE_CONTROL:
            self.negative_patients.append(sample)
        else:
            raise RuntimeError("Illegitimate MSI classification")

    def add_multiple_samples(self, samples: List[Sample]):
        for s in samples:
            self.add_patient(s)

    def __len__(self):
        return len(self.msi_patients)+len(self.mss_patients)+len(self.negative_patients)



class SamplesDB:
    def __init__(self):
        # self.db: ex. {'cancer_type_a': {MSS: [patient_a, patient_b,...]}
        self.db: Dict[str, SamplesSet] = dict()

    def cancer_types(self) -> List[str]:
        return list(self.db.keys())

    def add_patient(self, sample: Sample):
        if sample.cancer_type not in self.db:
            self.db[sample.cancer_type] = SamplesSet(sample.cancer_type, [], [], [])
        self.db[sample.cancer_type].add_patient(sample)

    def add_multiple_samples(self, samples: List[Sample]):
        for s in samples:
            self.add_patient(s)

    def get_all_patients(self) -> List[SamplesSet]:
        return [patient_set for patient_set in self.db.values()]

    def __len__(self):
        return sum([len(patient_set) for patient_set in self.db.values()])