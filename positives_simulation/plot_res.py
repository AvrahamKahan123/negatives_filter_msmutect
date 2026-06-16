import os
import numpy as np
import matplotlib.pyplot as plt
purities = [round(x, 2) for x in np.arange(0.05, 1.05, 0.05)]

d = "/home/avraham/MaruvkaLab/msmutect_postprocessing/data/positives_simple_simulation"

vmin, vmax = 0, 300
# for p in purities:
#     current_file = os.path.join(d, f"purity={p}.npy")
#     x = np.load(current_file)
#     x[x > vmax] = vmax
#
#     plt.figure()
#     im = plt.imshow(x, vmin=vmin, vmax=vmax, cmap='viridis')
#     name = f"purity={p}"
#     plt.title(name)
#     plt.colorbar(im, label="Intensity")  # ← adds colorbar
#     plt.ylabel("Reference_distribution")
#     plt.xlabel("Non_reference_distribution")
#     plt.savefig(output_path, dpi=300, bbox_inches='tight')
#     plt.close(fig)
#
#     print(f"Saved: {output_path}")
output_dir="."
for p in purities:
    current_file = os.path.join(d, f"purity={p}.npy")
    x = np.load(current_file)
    x[x > vmax] = vmax

    fig, ax = plt.subplots()
    im = ax.imshow(x, vmin=vmin, vmax=vmax, cmap='viridis')
    ax.set_title(f"purity={p}")
    cbar = fig.colorbar(im, ax=ax, label="Intensity")

    output_path = os.path.join(output_dir, f"purity_{p}.png")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)

    print(f"Saved: {output_path}")
