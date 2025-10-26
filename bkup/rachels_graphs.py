from typing import List, Tuple

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def bar_plot(df: pd.DataFrame):
    # Example data
    metrics = ["Intractability", "Centrality", "Irresolvability"]  # 3 groups
    israeli_sample = df[df["Site"]=="Israel"]
    us_sample = df[df["Site"]=="US"]
    intractable_column = 'Overall_Intractability7'
    central_column = 'Perceived_Centrality'
    irresolvable_column = 'Perceived_Irresolvability'
    values_israel = [israeli_sample[intractable_column].mean(), israeli_sample[central_column].mean(), israeli_sample[irresolvable_column].mean()]
    values_us = [us_sample[intractable_column].mean(), us_sample[central_column].mean(), us_sample[irresolvable_column].mean()]
    x = np.arange(len(metrics))  # [0, 1, 2]
    width = 0.35  # width of each bar

    fig, ax = plt.subplots(figsize=(7, 5))

    # Bars
    israel_bars = ax.bar(x - width / 2, values_israel, width, label="Israel")
    us_bars = ax.bar(x + width / 2, values_us, width, label="US", color="red")

    # Labels and title
    ax.set_ylabel("Means of Intractability Measures", fontsize=13, fontweight="bold")
    ax.set_yticks(list(range(0, 7)))
    ax.set_xticks(x)
    ax.set_xticklabels(metrics)
    ax.legend(fontsize=12)

    # Optional: add numbers on top of bars
    for bars in [israel_bars, us_bars]:
        for bar in bars:
            height = bar.get_height()
            ax.annotate(f"{height:.2f}",
                        xy=(bar.get_x() + bar.get_width() / 2, height),
                        xytext=(0, 3),  # 3 points above
                        textcoords="offset points",
                        fontsize=13,
                        fontweight="bold",
                        ha="center", va="bottom")
    ax.tick_params(axis='x', labelsize=13)
    ax.tick_params(axis='y', labelsize=13 )
    ax = plt.gca()
    for label in ax.get_xticklabels():
        label.set_fontweight("bold")
    # for label in ax.get_yticklabels():
    #     label.set_fontweight("bold")

    plt.tight_layout()
    plt.show()



def turn_df_to_points(df: pd.DataFrame, site: str) -> List[Tuple]:
    central_column = 'Perceived_Centrality'
    irresolvable_column = 'Perceived_Irresolvability'
    relevant_df = df[df["Site"]==site]
    x_coords = list(relevant_df[central_column])
    y_coords = list(relevant_df[irresolvable_column])
    return [(x, y) for x,y in zip(x_coords, y_coords)]


def plot_two_scatter(df: pd.DataFrame):
    """
    Plot two sets of points on a scatter plot.

    Parameters:
        points1 (list of tuple): list of (x, y) coordinates for first set
        points2 (list of tuple): list of (x, y) coordinates for second set
    """
    # Unpack coordinates
    central_column = 'Perceived_Centrality'
    irresolvable_column = 'Perceived_Irresolvability'
    israel_points = turn_df_to_points(df, "Israel")
    us_points = turn_df_to_points(df, "US")
    x1, y1 = zip(*israel_points)
    slope_israel, intercept_israel = np.polyfit(x1, y1, 1)
    x2, y2 = zip(*us_points)
    slope_us, intercept_us = np.polyfit(x2, y2, 1)
    plt.figure(figsize=(6, 6))
    # plt.scatter(x1, y1, color="blue", label="Israel")
    # plt.scatter(x2, y2, color="red", label="US")
    plt.scatter(x1, y1, color="blue", alpha=0.6, marker="o", label="Israel")
    plt.scatter(x2, y2, color="red", alpha=0.6, marker="x", label="US")
    x = np.arange(1, 8)
    plt.plot(x, slope_israel * x + intercept_israel, color="blue", linewidth=3)
    plt.plot(x, slope_us * x + intercept_us, color="red", linewidth=3 )

    plt.xlabel("Perceived Centrality", fontsize=14, fontweight="bold")
    plt.ylabel("Perceived Irresolvability", fontsize=14, fontweight="bold")
    plt.xticks(list(range(1, 8)), fontsize=12)
    plt.yticks(list(range(1, 8)), fontsize=12)
    # for label in plt.get_xticklabels()+plt.get_yticklabels():
    #     label.set_fontweight("bold")
    plt.legend(fontsize=16)
    plt.grid(True)
    plt.show()


df = pd.read_csv("rachels.csv")
print(df["Site"].value_counts())
print(df.columns)
plot_two_scatter(df)
