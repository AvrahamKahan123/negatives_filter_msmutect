from dataclasses import dataclass

from positives_simulation.Distribution import Distribution
from positives_simulation.RandomRequest import RandomRequest


@dataclass
class Response:
    request: RandomRequest
    succeeded : bool
    message : str = ""
    distribution: Distribution = None