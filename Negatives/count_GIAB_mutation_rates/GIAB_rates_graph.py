import math

import matplotlib.pyplot as plt


a={1: 8711, 2: 2428, 3: 108, 4: 236, 5: 57, 6: 9, 7: 5, 8: 0, 9: 0, 10: 1, 11: 0, 12: 0, 13: 0, 14: 0, 15: 0, 16: 0, 17: 0, 18: 0, 19: 0, 20: 0, 21: 0, 22: 0, 23: 0, 24: 0, 25: 0, 26: 0, 27: 0, 28: 0, 29: 0}

# counts in locus file
# b={1: 21849770, 2: 1709559, 3: 2283470, 4: 1061676, 5: 418856, 6: 194143, 7: 68969, 8: 41052, 9: 21042, 10: 17854, 11: 10476, 12: 9085, 13: 6528, 14: 6133, 15: 6794}
c=[1031530988,75409324,106593894,46645515,17652307,7702596,2295465,1141426,405644,335737,132851,148281,72459,74096,75027]
callable_bases = [6253343517,969485669,1214282503,776257381,339979931,168643542,59073891,35205493,14220536,13507271,5564049,6917782,3711353,3968225,4327529]
b=dict()
for i in range(len(c)):
    b[i+1]=callable_bases[i]
import matplotlib.pyplot as plt
import pandas as pd

# -----------------------------------------------------------------------------
# Data Preparation
# -----------------------------------------------------------------------------
# a = {
#     1: 8711,
#     2: 2428,
#     3: 108,
#     4: 236,
#     5: 57,
#     6: 9,
#     7: 5,
#     8: 0,
#     9: 0,
#     10: 1,
#     11: 0,
#     12: 0,
#     13: 0,
#     14: 0,
#     15: 0,
#     16: 0,
#     17: 0,
#     18: 0,
#     19: 0,
#     20: 0,
#     21: 0,
#     22: 0,
#     23: 0,
#     24: 0,
#     25: 0,
#     26: 0,
#     27: 0,
#     28: 0,
#     29: 0,
# }
#
# b = {
#     1: 21849770,
#     2: 1709559,
#     3: 2283470,
#     4: 1061676,
#     5: 418856,
#     6: 194143,
#     7: 68969,
#     8: 41052,
#     9: 21042,
#     10: 17854,
#     11: 10476,
#     12: 9085,
#     13: 6528,
#     14: 6133,
#     15: 6794,
# }

# Calculate False Positive Rate (FPR = a / b) for matching keys
motif_sizes = sorted(set(a.keys()).intersection(set(b.keys())))
fpr = [math.log10(a[k]/b[k]+1e-10) / b[k] for k in motif_sizes]

# -----------------------------------------------------------------------------
# Journal Figure Configuration (Publication Quality)
# -----------------------------------------------------------------------------
# Configure Matplotlib rcParams for journal standards
plt.rcParams["font.sans-serif"] = "Arial"
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["pdf.fonttype"] = 42  # Keep fonts as vector paths for EPS/PDF
plt.rcParams["ps.fonttype"] = 42

fig, ax = plt.subplots(figsize=(4.5, 3.5), dpi=300)

# Plot Data
ax.plot(
    motif_sizes,
    fpr,
    marker="o",
    color="#1F4E79",  # Classic journal navy blue
    linewidth=1.25,
    markersize=4.5,
    markerfacecolor="#1F4E79",
    markeredgecolor="none",
)

# -----------------------------------------------------------------------------
# Axes & Typography Styling
# -----------------------------------------------------------------------------
ax.set_xlabel("Motif Size", fontsize=10, fontweight="bold", labelpad=6)
ax.set_ylabel("FPR (Mutations/Callable Bases)", fontsize=10, fontweight="bold", labelpad=6)

ax.set_xticks(motif_sizes)
ax.set_xlim(0.5, 15.5)
ax.set_ylim(0, max(fpr) * 1.08)

# Configure ticks according to journal guidelines (inward ticks, clean frame)
ax.tick_params(
    axis="both",
    which="major",
    labelsize=8,
    direction="in",
    length=4,
    width=0.8,
    top=True,
    right=True,
)

# Set clean bounding box line width
for spine in ax.spines.values():
    spine.set_linewidth(0.8)
    spine.set_color("#000000")

plt.tight_layout()

# Save in high-res PNG and publication-ready vector format (PDF)
plt.savefig("fpr_vs_motif_size.png", dpi=300)
plt.savefig("fpr_vs_motif_size.pdf", format="pdf")
plt.show()