import sys


def ammend_line(line: str):
    fields = line.split("\t")
    return f"chr{fields[0]}:{fields[1]}-{fields[2].rstrip()}"


with open(sys.argv[1], 'r') as croc:
    l = croc.readlines()
regions = [ammend_line(line) for line in l]
with open("get_mutated_regions.sh", 'w+') as croc:
    croc.write("#!/bin/bash\n")
    croc.write(f"samtools view -H -b ${{1}} {' '.join(regions)} > ${{2}}")