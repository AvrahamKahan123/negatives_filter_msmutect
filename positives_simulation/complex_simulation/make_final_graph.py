import pandas as pd
import matplotlib.pyplot as plt


def plot_dict_lines(data):
    """
    data: dict with 5 keys, each value is a list of length 19
    """
    plt.figure(figsize=(10, 6))

    for key, values in data.items():
        plt.plot(values, marker='o', label=key)

    plt.xlabel("Index")
    plt.ylabel("Value")
    plt.title("Lines for Each Key")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def main():
    a={'TCGA-A6-5661': {10: 108, 15: 1735, 20: 11741, 25: 41638, 30: 98327, 35: 173006, 40: 251273, 45: 319894, 50: 373707, 55: 413944, 60: 442503, 65: 463312, 70: 478156, 75: 489295, 80: 497325, 85: 503205, 90: 507344, 95: 510207, 100: 512109}, 'TCGA-AJ-A3BH': {10: 115, 15: 1413, 20: 8432, 25: 30220, 30: 72177, 35: 129533, 40: 190930, 45: 247423, 50: 293540, 55: 330052, 60: 357797, 65: 378497, 70: 394080, 75: 406220, 80: 414679, 85: 421032, 90: 425493, 95: 428307, 100: 430255}, 'TCGA-AP-A05N': {10: 124, 15: 1445, 20: 7926, 25: 24171, 30: 52098, 35: 86903, 40: 123528, 45: 156662, 50: 183910, 55: 204617, 60: 220969, 65: 232589, 70: 241567, 75: 248765, 80: 253867, 85: 258217, 90: 261406, 95: 263780, 100: 265581}, 'TCGA-FI-A2D4': {10: 101, 15: 1378, 20: 8121, 25: 26758, 30: 60704, 35: 105895, 40: 153922, 45: 198693, 50: 235769, 55: 265107, 60: 287940, 65: 304662, 70: 317201, 75: 327176, 80: 334218, 85: 339862, 90: 343832, 95: 346853, 100: 349066}, 'TCGA-OR-A5LB': {10: 37, 15: 446, 20: 2568, 25: 8637, 30: 19560, 35: 34320, 40: 50942, 45: 66620, 50: 80359, 55: 92002, 60: 101493, 65: 108972, 70: 115294, 75: 120002, 80: 123862, 85: 127007, 90: 129350, 95: 130875, 100: 132221}}
    total_mut_counts = {'TCGA-A6-5661': 798241, 'TCGA-AJ-A3BH': 667610, 'TCGA-AP-A05N': 437981, 'TCGA-FI-A2D4': 561264,
                        'TCGA-OR-A5LB': 242120 }
    case_lines = {}
    for case in a.keys():
        l = []
        for i in range(10, 105, 5):
            l.append(a[case][i]/total_mut_counts[case])
        case_lines[case] = l
    print(case_lines)
    plot_dict_lines(case_lines)
if __name__ == "__main__":
    main()
# v=pd.DataFrame(a)
# print(v)
