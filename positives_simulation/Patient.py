import os
from typing import List, Tuple

from positives_simulation.DBwriteRequest import DBwriteRequest
from positives_simulation.Distribution import Distribution
from positives_simulation.PatientDBengine import PatientDBengine


class Patient:
    def __init__(self, patient_name: str, patient_path: str, create=False) -> None:
        self.patient_name = patient_name
        self.patient_path = patient_path
        self.homozygous_reference_dir = os.path.join(self.patient_path, "Homozygous_Reference")
        self.heterozygous_reference_dir = os.path.join(self.patient_path, "Heterozygous_Reference")
        self.non_reference_homozygous_dir = os.path.join(self.patient_path, "Non_Reference_Homozygous")
        if create:
            for folder in [self.homozygous_reference_dir, self.heterozygous_reference_dir, self.non_reference_homozygous_dir]:
                os.makedirs(folder, exist_ok=True)
        self.homozygous_reference_db = PatientDBengine(self.homozygous_reference_dir, create=create)
        self.heterozygous_reference_db = PatientDBengine(self.heterozygous_reference_dir, create=create)
        self.non_reference_db = PatientDBengine(self.non_reference_homozygous_dir, create=create)

    def get_random_homozygous_reference(self, reference_size: int) -> Distribution:
        return self.homozygous_reference_db.get_random_entry([reference_size])

    def get_random_heteroyzgous_reference(self, reference_size: int) -> Distribution:
        return self.heterozygous_reference_db.get_random_entry([reference_size])

    def get_random_non_reference(self, reference_size: int, non_reference_size: int) -> Distribution:
        return self.non_reference_db.get_random_entry([reference_size, non_reference_size])

    def size(self):
        raise self.homozygous_reference_db.size() + self.heterozygous_reference_db.size() + self.non_reference_db.size()

    def divide_write_requests(self, write_requests: List[DBwriteRequest]) -> Tuple[List[DBwriteRequest], List[DBwriteRequest], List[DBwriteRequest], List[DBwriteRequest]]:
        (homozygous_reference_requests, heterozygous_reference_requests,
         non_reference_requests, heterozygous_non_reference_requests) = [], [], [], []
        for write_request in write_requests:
            if write_request.homozygous:
                if write_request.reference:
                    homozygous_reference_requests.append(write_request)
                else:
                    non_reference_requests.append(write_request)
            else:
                if write_request.reference:
                    heterozygous_reference_requests.append(write_request)
                else:
                    heterozygous_non_reference_requests.append(write_request)
        return homozygous_reference_requests, heterozygous_reference_requests, non_reference_requests, heterozygous_non_reference_requests

    def add_distributions(self, write_requests: List[DBwriteRequest]) -> None:
        homozygous_reference_requests, heterozygous_reference_requests, non_reference_requests, _ = self.divide_write_requests(write_requests)
        self.homozygous_reference_db.write_entries(homozygous_reference_requests)
        self.heterozygous_reference_db.write_entries(heterozygous_reference_requests)
        self.non_reference_db.write_entries(non_reference_requests)

