
import csv
from itertools import islice

a={i:0 for i in range(5, 41)}
print(a)
with open("/home/avraham/MaruvkaLab/msmutect_development/data/GRCh38.d1.vd1_1to15_repetitive_loci_sorted_fixed",
          "r", newline="") as f:
    reader = csv.DictReader(f, delimiter="\t", fieldnames=["chrom", "rep_length"]+[f"col{x}" for x in range(4)]+["NUM_REPEATS"]+[f"col{i+7}" for i in range(8)])
    # for row in islice(reader, 50):
    for row in reader:
        # Each row is a dict: {column_name: value}
        repeat_length = int(float(row["rep_length"]))
        if repeat_length==1:
            num_repeats = int(float(row["NUM_REPEATS"]))
            if num_repeats<41:
                a[num_repeats]+=1
print(a)