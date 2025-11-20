#!/usr/bin/env -S python 
import sys, os
sys.path.append('../py')
sys.stdout.reconfigure(line_buffering=True)

from netcdf import read_netcdf, extract_dataframe
import numpy as np
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt
# matplotlib.use('Agg')  ## Noninteractive

fn = f"{os.path.dirname(__file__)}/my.mplstyle"
plt.style.use(fn)


## Dataset manipulation
src_cdf = "run/cache/ds_point.cd" 
ds = read_netcdf(src_cdf)

df = extract_dataframe( ds, "UY" )
df["strain"] = - df.UY / 10


# Reference data
ref_df = pd.read_csv("raw_data/poiate-2012-transient.csv", sep="\t")
ref_df["time_in_days"] = ref_df.time_h/24
# Offset to avoid transients and lab issues in early times
ref_df["strain"] -= 0.014

## Plot simulation
fig, ax = plt.subplots( figsize=(4,2.5))
_df = df[::50]
ax.plot(_df.time_in_days, _df.strain, c='k', lw=2, label = "Model")

# Plot raw data - interpolate for a better plot
from scipy.interpolate import interp1d
f_ = interp1d( ref_df.time_in_days , ref_df.strain, fill_value='extrapolate')
t_ = np.linspace(-10,90,30)
s_ = f_(t_)
ax.scatter(t_, s_, c='r', marker='x', lw=1.5, label= "Experimental")

ax.set_xlim( 0, 90 )
ax.set_ylim( 0, 0.2 )
ax.legend( bbox_to_anchor=(0.02, 0.98), loc="upper left")

ax.set_xlabel('Time (days)')
ax.set_ylabel('Strain')

import plot_util
plot_util.savefig( fig, "png/transient.png" )


# plt.savefig('uy.png', dpi=300, bbox_inches='tight')
# plt.close()

