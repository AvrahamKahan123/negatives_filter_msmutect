import os
from typing import Dict, List

from positives_simulation.DBwriteRequest import DBwriteRequest
from positives_simulation.Distribution import Distribution
from positives_simulation.Patient import Patient
from positives_simulation.RandomRequest import RandomRequest
from positives_simulation.Response import Response


class RandomAllelesDB:
    # AllelesDB doesn't care about ACID, or nothing
    def __init__(self, db_path: str):
        # db_path is a director
        self.db_path = db_path
        if not os.path.exists(self.db_path):
            os.makedirs(self.db_path)
        self.patients = os.listdir(self.db_path)

        self._loaded_patients: Dict[str, Patient] = dict()


    def loaded_patients(self):
        return list(self._loaded_patients.keys())

    def load_patient_files(self, patient_name: str, patient_path: str):
        self._loaded_patients[patient_name] = Patient(patient_name, patient_path)
    
    def load_patient(self, patient_name: str):
        patient_path = os.path.join(self.db_path, patient_name)
        if patient_name not in self.patients:
            raise RuntimeError(f"patients {patient_name} does not exist")
        self.load_patient_files(patient_name, patient_path)

    def unload_patient(self, patient_name: str):
        if patient_name not in self.loaded_patients:
            raise RuntimeError(f"Patient {patient_name} is not currently loaded")
        del self._loaded_patients[patient_name]

        
    def memory_usage(self) -> int:
        return sum([patient.size() for patient in self._loaded_patients.values()])

    def get(self, request: RandomRequest) -> Response:
        if request.patient_name not in self._loaded_patients:
            return Response(succeeded=False, message="Patient {request.patient_name} is not currently loaded",
                            request=request)

        patient = self._loaded_patients[request.patient_name]
        return_distribution = self.query_patient(patient, request)
        if len(return_distribution.repeat_lengths) == 0:
            return Response(succeeded=False, message="No distributions matched request", request=request)
        else:
            return Response(succeeded=True, distribution=return_distribution, request=request)

    def query_patient(self, patient: Patient, request: RandomRequest) -> Distribution:
        if request.follows_reference():
            if request.homozygous:
                return patient.get_random_homozygous_reference(request.reference_size)
            else: # heterozygous
                return patient.get_random_heteroyzgous_reference(request.reference_size)
        else: # doesn't follow reference
            return patient.get_random_non_reference(request.reference_size, request.true_size)

    def add_patient(self, patient_name: str, write_requests: List[DBwriteRequest]) -> None:
        patient_path = os.path.join(self.db_path, patient_name)
        if not os.path.exists(patient_path):
            # raise RuntimeError(f"Patient {patient_name} already exists")
            os.makedirs(patient_path)
        patient = Patient(patient_name, patient_path, create=True)
        patient.add_distributions(write_requests)

        