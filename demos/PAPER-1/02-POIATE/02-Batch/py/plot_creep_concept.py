#!/usr/bin/env -S python 

import sys, os
sys.path.append('../py')

import numpy as np
import matplotlib.pyplot as plt

from plot_util import savefig
cwd = os.path.dirname(__file__)
fn = f"{cwd}/my.mplstyle"
plt.style.use(fn)


# Define stress range (MPa)
stress = np.logspace(0, 1.8, 500)  # 0.1 to 100 MPa


# sig = sig_MPa * 1e6
T = 50 + 273
fud = 1.08  ##  A Fudge factor (see arrhenius)

Q = 51600
R = 8.314

## Steady State
sig0ss = [ 10e6 ,   50e6  ]
N    = [    8.6 ,   .7    ]

arr = np.exp( -(Q/R/T)**fud  )

# Steady State
# eps_ss_rate = 0
# for i in range(2) : 
#     eps_ss_rate += np.exp( -(Q/R/T)**fud  ) * ( sig/sig0ss[i] )**N[i]


epsilon_dot_ps = arr*(stress/50) 

n = 7.0  
epsilon_dot_dc = arr*(stress/10)**n

# Combined creep (sum of mechanisms)
epsilon_dot_total = epsilon_dot_ps + epsilon_dot_dc

fig, ax = plt.subplots(figsize=(4, 2.5))

ax.loglog(stress, epsilon_dot_total, 'k-', linewidth=1.0, zorder=10)

# Add asymptotes as dashed lines
# Low stress asymptote (pressure solution dominates)
ax.loglog(stress, epsilon_dot_ps, 'r--', linewidth=0.7, alpha=0.7,
          label='Pressure solution asymptote')

ax.loglog(stress, epsilon_dot_dc, 'b--', linewidth=0.7, alpha=0.7,
          label='Dislocation creep asymptote')

# Labels
ax.set_xlabel('Deviatoric Stress (MPa)')
ax.set_ylabel('Strain rate (1/s)')
from matplotlib.ticker import ScalarFormatter
ax.xaxis.set_major_formatter(ScalarFormatter())

# Legend - only for asymptotes
ax.legend(loc='upper left', bbox_to_anchor=(0.02, 0.98))
ax.set_xticks([1,4,10,40])

# Set axis limits
ax.set_xlim([1,60])
ax.set_ylim([1e-14, 1e-5])

plt.tight_layout()
savefig(fig,'png/salt_creep_mechanisms.png')
print("Plot saved successfully")

plt.show()

