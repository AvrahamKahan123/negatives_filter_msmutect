from dataclasses import dataclass

@dataclass
class RandomRequest:
    def __init__(self, patient_name: str, homozygous: bool, reference_size: int, true_size: int = None):
        self.patient_name = patient_name
        self.homozygous = homozygous
        self.reference_size = reference_size
        if true_size is None:
            self.true_size = reference_size
        else:
            self.true_size = true_size

        if (not self.homozygous) and (not self.follows_reference()):
            raise RuntimeError(f"Invalid Request: Heterozygous, non-reference alleles are not allowed")

    def follows_reference(self):
        return self.true_size == self.reference_size