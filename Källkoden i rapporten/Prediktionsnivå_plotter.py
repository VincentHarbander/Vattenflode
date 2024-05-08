# målet är att plotta prediktion o projektion i det här programmet
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
from collections import namedtuple

import matplotlib.colors as mcolors

# Import data
river_directory = r"R/CSVs/bm0"
River = namedtuple("River", ["name", "data"])
rivers = []
# Itererar över alla filer i nedladdningsmappen
for filename in os.listdir(river_directory):
    filepath = os.path.join(river_directory, filename)
    # Ifall filer ändras / raderas under tiden som programmet körs.
    if not os.path.isfile(filepath):
        continue

    river_name = filename.split(".")[0]
    rivers.append(
        River(
            name=river_name,
            data=pd.read_csv(
                filepath,
                delimiter=",",
                decimal=".",
            ),
        )
    )


# för nu använder vi ".1" för low och ingenting för high.
low_column_char = ""
high_column_char = ".1"
low_column_names = [
    str(round(i, 2)) + low_column_char for i in list(np.arange(0.01, 1, 0.01))
]
high_column_names = [
    str(round(i, 2)) + high_column_char for i in list(np.arange(0.01, 1, 0.01))
]
column_names = low_column_names + high_column_names


for river in rivers:
    cmap = plt.get_cmap("plasma")
    norm = mcolors.Normalize(vmin=0, vmax=len(column_names))

    river.data["Unnamed: 0"] = river.data["Unnamed: 0"] + 2024  # gör om till år
    for i in range(len(high_column_names) - 1):
        # Adjust alpha value to distribute transparency equally among lines
        alpha = 1.0 / len(column_names)
        color = cmap(norm(i))  # Get color from colormap
        plt.fill_between(
            river.data["Unnamed: 0"],
            river.data[low_column_names[i]],
            river.data[low_column_names[i + 1]],
            # alpha=alpha,
            color=color,
        )
        plt.fill_between(
            river.data["Unnamed: 0"],
            river.data[high_column_names[i]],
            river.data[high_column_names[i + 1]],
            # alpha=alpha,
            color=color,
        )

    plot_lines = True
    if plot_lines:
        # Lägger till linjer för 95% konfidensintervall
        plt.plot(
            river.data["Unnamed: 0"],
            river.data["0.95"],
            color="black",
            linestyle="dashed",
            linewidth=1,
        )
        plt.plot(
            river.data["Unnamed: 0"],
            river.data["0.95.1"],
            color="black",
            linestyle="dashed",
            linewidth=1,
        )
        # och för 99% konfidensintervall
        plt.plot(
            river.data["Unnamed: 0"],
            river.data["0.99"],
            color="black",
            linestyle="dashed",
            linewidth=1,
        )
        plt.plot(
            river.data["Unnamed: 0"],
            river.data["0.99.1"],
            color="black",
            linestyle="dashed",
            linewidth=1,
        )
        # och för medelvärdet
        plt.plot(
            river.data["Unnamed: 0"],
            river.data["mean"],
            color="black",
            linewidth=1,
        )

    # TODO: colormap.
    # sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    # sm.set_array([])
    # plt.colorbar(sm, label="Colorbar Title")
    plt.ylabel("Prediktion av vattenföring (m^3/s)")
    plt.xlabel("År")
    savepath = r"Python/Grafresultat/Prediktionsplottar/"
    savepath = os.path.join(savepath, f"{river.name}.png")
    plt.savefig(savepath, bbox_inches="tight")
    plt.close()
    # plt.show()
    # raise NotImplementedError  # kör det här när vi faktiskt sparar grejerna
