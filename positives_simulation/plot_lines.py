# import pandas as pd
# import matplotlib.pyplot as plt
#
# # Load the CSV
# df = pd.read_csv("../results/complex_simulation_stats/percentages_table.csv")
#
# # Plot all columns
# plt.figure(figsize=(10, 6))
# for column in df.columns:
#     plt.plot(df.index, df[column], label=column)
#
# plt.xlabel("Index")
# plt.ylabel("Value")
# plt.title("All Columns as Lines")
# plt.legend()
# plt.grid(True)
# plt.tight_layout()
# plt.show()


import pandas as pd
import matplotlib.pyplot as plt

# Load CSV
df = pd.read_csv("../results/complex_simulation_stats/percentages_table.csv")

# Use the first column as x, and the rest as y
x = df.iloc[:, 0]
y_cols = df.columns[1:]

# Plot each column as a separate line
plt.figure(figsize=(10, 6))
for col in y_cols:
    plt.plot(x, df[col], label=col)

plt.xlabel(df.columns[0])  # label x-axis with the first column name
plt.ylabel("Mutation Called Percentage")
plt.title("Mutation Percentages Per Different Purities")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
