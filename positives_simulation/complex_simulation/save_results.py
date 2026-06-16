
def extract_segment(line):
    # Find last underscore
    pos1 = line.rfind("_")
    if pos1 == -1:
        return None

    # Find first colon AFTER the last underscore
    pos2 = line.find(":", pos1)
    if pos2 == -1:
        return None

    return line[pos1+1 : pos2]

def make_dict_for_case():
    return {5*i: -1 for i in range(2,21)}

def main():
    with open("results.txt", 'r') as f:
        lines = f.readlines()
    cases = ["TCGA-A6-5661", "TCGA-AJ-A3BH", "TCGA-AJ-A3BH", "TCGA-AP-A05N", "TCGA-FI-A2D4", "TCGA-OR-A5LB"]
    results = {case: make_dict_for_case() for case in cases}
    print(results)
    for line in lines:
        case = line[:12]
        purity = float(extract_segment(line))
        int_purity = int(purity*100)
        num_mutations = int(line[line.rfind(":")+1:].strip())
        results[case][int_purity] = num_mutations
    print(results)

if __name__ == "__main__":
    main()
