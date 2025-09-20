#!/usr/bin/env -S python -i
import matplotlib.colors as mcolors
import matplotlib.cm as cm

import xarray as xr
ds = xr.open_dataset("run/cdf/plane_0.cd")

lab = ds.Coord.attrs['components'].split() 
ds = ds.assign_coords(vec3_comp=lab)


# lab = ds.Stress.attrs['components'].split() 
# ds = ds.assign_coords(ten9_comp=lab)
# print(ds.Stress.sel(ten9_comp='xy'))

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import xarray as xr
import numpy as np

# Suppose your dataset is already loaded as ds

# Extract Y and Z coordinates (fixed in time)
y = ds.Coord[:, ds.vec3_comp.values.tolist().index("y")].values
z = ds.Coord[:, ds.vec3_comp.values.tolist().index("z")].values

# Extract pressure (time, point_idx)
D = { 'vname':"Delta_T" , 'min':-10, 'max':0 }
D = { 'vname':"Delta_P" , 'min':None, 'max':1e7 }
VARNAME = D['vname']
vals = ds[VARNAME].values 
times = ds.time.values

fig, ax = plt.subplots(figsize=(6, 5))
sc = ax.scatter(y, z, c=vals[0], s=20, cmap="coolwarm",  vmin=D['min'], vmax=D['max'])
cb = plt.colorbar(sc, ax=ax, label=VARNAME)
time_text = ax.text(0.02, 0.95, "", transform=ax.transAxes)

ax.set_xlabel("Y")
ax.set_ylabel("Z")
ax.set_title(f"{VARNAME}")
ax.invert_yaxis()

def update(frame):
    sc.set_array(vals[frame])
    time_text.set_text(f"t = {times[frame]:.2f}")
    return sc, time_text

ani = FuncAnimation(fig, update, frames=len(times), interval=200, blit=False)

plt.show()

