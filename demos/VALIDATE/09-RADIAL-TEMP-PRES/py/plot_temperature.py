#!/usr/bin/env -S python

import pandas as pd
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt

# Load CSV file (assumes it has no header; adjust if needed)
df = pd.read_csv('t.csv', header=None)

# Extract columns (adjust index if header exists)
x = df[1]  # Second column
y = df[2]  # Third column
z = df[3]  # Fourth column, for color

contour = plt.scatter(x, y, c=z, s=10, linestyle="None")

plt.colorbar(contour, label='Interpolated Value (Column 4)')
plt.xlabel('X (Column 2)')
plt.ylabel('Y (Column 3)')
plt.title('Interpolated Colormap from CSV Data')

plt.xlim( 0, 1 )
plt.tight_layout()
plt.show()
