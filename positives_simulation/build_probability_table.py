import numpy as np

def row_cumulative_sum(matrix):
    result = []
    for row in matrix:
        cum = 0
        new_row = []
        for x in row:
            cum += x
            new_row.append(cum)
        result.append(new_row)
    return result

collated_npy = "/home/avraham/MaruvkaLab/msmutect_postprocessing/positives_simulation/simple_simulation/collated.npy"
x=np.load(collated_npy)
sum_per_column =(x.sum(axis=1).clip(min=1))
xs = x/sum_per_column[:, None]#np.divide(x, sum_per_column, axis=1)
cumsum = np.array(row_cumulative_sum(xs))
print(cumsum)
print(cumsum[4])
np.save("p_mut_map.npy", cumsum)
# print(x.shape)
# print(x.sum(axis=1))
# for i in range(len(x)):
#     current_row = x[i]
#     s = max(current_row.sum(), 1)
#     current_row_duplicated = current_row.copy()
#     normalized = current_row_duplicated/s
#     print(normalized)