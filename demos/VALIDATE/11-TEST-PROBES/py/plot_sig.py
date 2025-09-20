#!/usr/bin/env -S python -i

import pandas as pd
import matplotlib.pyplot as plt

# Filter for timestep=30 and variable=STOTXX
# TS = 23
# df = pd.read_csv("run.0/csv/plane_0.csv", sep="\t")
# df = df[(df["Timestep"] == TS)]
# df = df[(df["Var"] == "invarQ") | (df.Var == "invarPeff")]
# print("Pivoting ...")
# df = df.pivot_table(index=["Timestep", "X", "Y", "Z"], columns="Var", values="Value").reset_index()
# print("Writing pickle ...")
# df.to_pickle("run.0/pivot.pkl")

df = pd.read_pickle("run.0/pivot.pkl")

# # Make the scatter plot
plt.scatter(df.invarPeff, df.invarQ )
plt.xlabel("P'")
plt.ylabel("Q")
plt.gca().invert_yaxis()

# plt.show()
