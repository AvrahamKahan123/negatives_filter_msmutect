import csv
import sys
from collections import Counter

with open(sys.argv[1], newline="") as f:
    counts = Counter(len(row["PATTERN"].strip())
                     for row in csv.DictReader(f, delimiter="\t"))

for length in sorted(counts):
    print(f"{length}\t{counts[length]}")