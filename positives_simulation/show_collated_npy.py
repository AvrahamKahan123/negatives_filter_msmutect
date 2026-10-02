import numpy as np
import matplotlib.pyplot as plt

collated_npy_fp = "/home/avraham/MaruvkaLab/msmutect_postprocessing/positives_simulation/simple_simulation/collated.npy"
x = np.load(open(collated_npy_fp, "rb"))


## create a graph of the diagonal line one over the main diagonal minus the one
## under it, taken row by row, normalized by its sum.
## entry i is C[i, i+1] - C[i, i-1], so the first entry is C[1,2] - C[1,0].
## rows 0 and 40 are dropped: they have no lower / no upper neighbour.
n = x.shape[0]
rows = np.arange(1, n - 1)
above = x[rows, rows + 1]   # one repeat unit longer than the reference
below = x[rows, rows - 1]   # one repeat unit shorter than the reference

diff = above - below
normalized = diff / diff.sum()

plt.figure(figsize=(10, 5))
plt.plot(rows, normalized, marker="o", markersize=4)
plt.axhline(0, color="black", linewidth=0.8)
plt.xlabel("Reference_distribution")
plt.ylabel("(above - below) / sum")
plt.title("Expansion vs Contraction Bias, One Repeat Unit")
plt.xticks(np.arange(0, n, 2))
plt.grid(alpha=0.3)
plt.savefig("collated_offdiagonal_bias.png", dpi=300, bbox_inches="tight")
plt.show()
