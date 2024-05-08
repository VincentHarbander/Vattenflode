# Vincent
# 2024

# pdflatex main.tex --shell-escape

import numpy as np, pandas as pd

all_data = "R/All data.csv"

df = pd.read_csv(all_data, sep=";", encoding="utf-8")

result = ""
example = ""

with open("input.txt", mode="r", encoding="utf-8") as file:
    example = file.read()

for i in range(len(df)):
    curr = df.iloc[i,:]

    raw = curr["name"][:-4]  # Remove the .csv ending
    split_list = raw.split("_")
    underscoredrealname = "_".join(split_list[:-1])
    name = " ".join(split_list[:-1])
    station = split_list[-1]

    chunk = example.replace("underscored", raw).replace("realname", name).replace("stationnr", station)

    # These are replaced first because "lambda0" is a substring of "lambda0_se", which is inconvenient!
    chunk = chunk.replace(" lambda0_se", " \\num{" + str(curr["lambda0_se"]) + "}")
    chunk = chunk.replace(" lambda1_se", " \\num{" + str(curr["lambda1_se"]) + "}")
    chunk = chunk.replace(" lambda1_p", " \\num{" + str(curr["lambda1_p"]) + "}")
    chunk = chunk.replace(" lambda1_bh", " \\num{" + str(curr["lambda1_bh"]) + "}")

    for param in curr.keys():
        if param == "post_datapts" or str(curr[param]) == "nan":
            chunk = chunk.replace(" " + param, " " + str(curr[param]))
        elif param in ["yellow_level", "orange_level", "red_level"]:
            chunk = chunk.replace(param, "\\qty{" + str(curr[param]) + "}{\\meter\\cubed\\per\\second}")
        else:
            chunk = chunk.replace(" " + param, " \\num{" + str(curr[param]) + "}")

    result += chunk

with open("LaTeX.txt", mode="w", encoding="utf-8") as file:
    file.write(result)
