#!/usr/bin/env -S python -i

import pandas as pd
import numpy as np

# Load CSV
df = pd.read_csv("run/csv/plane_0.csv", sep="\t")

# # Pivot so that each variable becomes a column
pivot = df.pivot_table(index=["Timestep", "X", "Y", "Z"],
                       columns="Var", values="Value").reset_index()

# # Function to compute invariants for one timestep
# def compute_invariants(group):
#     # Effective stresses: σ' = σ_total - P
#     sxx = group["STOTXX"] - group["P"]
#     syy = group["STOTYY"] - group["P"]
#     szz = group["STOTZZ"] - group["P"]
#     sxy = group["STOTXY"]
#     sxz = group["STOTXZ"]
#     syz = group["STOTYZ"]

#     # Mean effective stress
#     p_eff = (sxx + syy + szz) / 3

#     # Deviatoric stresses
#     sxx_dev = sxx - p_eff
#     syy_dev = syy - p_eff
#     szz_dev = szz - p_eff

#     # J2 invariant
#     j2 = 0.5 * (
#         sxx_dev**2 + syy_dev**2 + szz_dev**2 +
#         2 * (sxy**2 + sxz**2 + syz**2)
#     )

#     q = np.sqrt(3 * j2)

#     return pd.DataFrame({
#         "p_eff": p_eff,
#         "q": q,
#         "q/p_eff": q / p_eff
#     })

# # Apply to each timestep
# results = pivot.groupby("Timestep").apply(compute_invariants).reset_index()

# # If you want one value per timestep (mean over points)
# summary = results.groupby("Timestep")[["p_eff", "q", "q/p_eff"]].mean()

# print(summary)
