import csv
from typing import Dict

import numpy as np

def purities():
    return [round(x, 2) for x in np.arange(0, 1.05, 0.05)]

class Case:
    def __init__(self, case_name: str, num_mutations: int):
        self.case_name = case_name
        self.num_mutations = num_mutations
        self.percentages_dict = {p: 0 for p in purities()}

    def add_percentage(self, purity: float, num_mutations: int):
        if purity not in self.percentages_dict:
            raise RuntimeError("Purity not found")
        self.percentages_dict[purity] = num_mutations/self.num_mutations



def main():
    TCGA_A6_5661 =  Case("TCGA-A6-5661", 830832)
    TCGA_AJ_A3BH = Case("TCGA-AJ-A3BH", 686049)
    TCGA_AP_A05N = Case("TCGA-AP-A05N", 469233)
    TCGA_FI_A2D4 = Case("TCGA-FI-A2D4", 592553)
    TCGA_OR_A5LB = Case("TCGA-OR-A5LB", 253848)
    cases: Dict[str, Case] = {case.case_name: case for case in [TCGA_A6_5661, TCGA_AJ_A3BH, TCGA_AP_A05N, TCGA_FI_A2D4, TCGA_OR_A5LB]}

    with open("../../results/complex_simulation_stats/total_mut_count_table.txt", "r") as file:

        reader = csv.DictReader(file, delimiter='\t')
        for row in reader:
            filename = row["Case"]
            case_name = filename[:12]
            purity = filename[filename.rfind("_")+1:]
            case = cases[case_name]
            case.add_percentage(float(purity), int(row["Total_mut_count"]))
    different_purities = purities()
    with open("../../results/complex_simulation_stats/percentages_table.csv", "w") as file:
        file.write("purity,"+",".join([case.case_name for case in cases.values()])+"\n")
        for p in different_purities:
            file.write(",".join([str(p)]+[str(case.percentages_dict[p]) for case in cases.values()])+"\n")


if __name__ == "__main__":
    main()
    


