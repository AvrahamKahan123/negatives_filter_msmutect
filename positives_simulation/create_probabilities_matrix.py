import numpy as np


def main():
    np_path = "/home/avraham/MaruvkaLab/msmutect_postprocessing/population_level_analysis/final_pop.npy"
    full_np = np.load(np_path)
    ret_np = np.zeros((41, 41))
    for i in range(len(full_np)):
        current_row = full_np[i]
        reference = current_row.argmax()
        current_row[reference] = 0
        row_num_mutated_samples = current_row.sum()
        if row_num_mutated_samples == 0:
            continue
        nonzero_idxs = current_row.nonzero()[0]
        nonzero_vals = current_row[nonzero_idxs]
        for idx, val in zip(nonzero_idxs, nonzero_vals):
            ret_np[reference, idx] += val
    np.save("simple_simulation/collated.npy", ret_np)


if __name__ == "__main__":
    main()