from dataclasses import dataclass
from typing import List

from positives_simulation.Distribution import Distribution


@dataclass
class DBwriteRequest:
    keys: List[int]
    homozygous: bool
    reference: bool
    distribution: Distribution