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

# Create a grid to interpolate onto
xi = np.linspace(min(x), max(x), 300)
yi = np.linspace(min(y), max(y), 300)
xi, yi = np.meshgrid(xi, yi)

# Interpolate z values onto the grid
zi = griddata((x, y), z, (xi, yi), method='cubic')  # use 'linear' or 'nearest' as alternatives

# Plot
plt.figure(figsize=(8, 6))
contour = plt.imshow(zi, extent=[x.min(), x.max(), y.min(), y.max()],
                     origin='lower', cmap='viridis', aspect='auto')

plt.colorbar(contour, label='Interpolated Value (Column 4)')
plt.xlabel('X (Column 2)')
plt.ylabel('Y (Column 3)')
plt.title('Interpolated Colormap from CSV Data')
plt.tight_layout()
plt.show()
