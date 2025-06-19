#!/usr/bin/env -S python

import os, time
import re
import numpy as np
import pandas as pd

# Parse the dat to capture the radial grid
dat_fn = "STARS/test3-explicit_well.dat"
temp_fn = "STARS/test3-refine-dates Temperature.txt"

# Precompile regex patterns for performance
patterns = {
    "GRID":    re.compile(r"^\s*\*?GRID\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+([0-9.]+)$"),
    "DI":      re.compile(r"^\s*\*?DI\s+(\S+)\s+(.+)$"),
    "DJ":      re.compile(r"^\s*\*?DJ\s+(\S+)\s+(.+)$"),
    "DK":      re.compile(r"^\s*\*?DK\s+(\S+)\s+(.+)$"),
    "DEPTH":   re.compile(r"^\s*\*?DEPTH\s+(\d+)\s+(\d+)\s+(\d+)\s+([0-9.]+)\s*$"),
}

data = { "DI" : [] , "DJ" : [], "DK" : [] }

state = "START"
prop = None
reading_mode = ""
n_grid = 0
idx = 0
depth = 0

#
#
#
def _append( l, vstr ) :
    l.extend( vstr.split() )

#
#
#
def _to_list( vstr ) :
    ret = []
    l1 = vstr.split()
    for s in l1 :
        r = re.compile( r"^(?:(\d+)?\*)?([0-9.]+)$" )
        m = r.match( s )
        if not m : 
            print("Error! => " + s)
            exit(-1)

        t = m.group(1)
        if not t : t = 1
        else : t = int(t)

        for i in range(t) :
            ret.append( float(m.group(2)) )

    return ret
        
#
#
#
with open(dat_fn, 'r') as file:
    ln = 0
    for line in file:
        line = line.strip()
        line = re.sub(r'\s*\*\*.*', '', line)
        ln += 1

        # An empty line resets the machine
        if not line :
            state = "START"
            prop = None
            continue

        if m := patterns["DEPTH"].match(line):
            depth = float( m.group(4) )
            state = "START"
            print(f"Depth:{depth}")
            continue

        if m := patterns["GRID"].match(line):
            t = m.group(1).lstrip("*").lower()
            if t != "radial" : 
                print("We only support radial grids!")
                exit
            ni = int(m.group(2))
            nj = int(m.group(3))
            nk = int(m.group(4))
            rwstr = m.group(5).lstrip("*")
            rw = float(m.group(6))
            state = "START"
#             print(f"ni:{ni} nj:{nj} nk:{nk} rwstr:{rwstr} rw:{rw}")
            continue

        if m := patterns["DI"].match(line):
            reading_mode = m.group(1).lstrip("*").lower()
            line = m.group(2)
            state = "READING"
            data["DI"] = np.zeros( (ni, nj, nk) )
            prop = data["DI"]
            idx=0

        if m := patterns["DJ"].match(line):
            reading_mode = m.group(1).lstrip("*").lower()
            line = m.group(2)
            state = "READING"
            data["DJ"] = np.zeros( (ni, nj, nk) )
            prop = data["DJ"]
            idx=0

        if m := patterns["DK"].match(line):
            reading_mode = m.group(1).lstrip("*").lower()
            line = m.group(2)
            state = "READING"
            data["DK"] = np.zeros( (ni, nj, nk) )
            prop = data["DK"]
            idx=0

        if state == "READING" :
            if reading_mode == "con" :
                v = float( line )
                prop[:,:,:] = v 

            vl = _to_list( line )
            for v in vl :
                if reading_mode == "ivar" : prop[idx,:,:] = v
                if reading_mode == "jvar" : prop[:,idx,:] = v
                if reading_mode == "kvar" : prop[:,:,idx] = v

                idx += 1

#
#  Calculate the cell coordinates
#

# Separate components
di = data["DI"]
dj = data["DJ"]
dk = data["DK"]

# Initialize coordinate arrays for cell centers
x = np.zeros((ni, nj, nk))
y = np.zeros((ni, nj, nk))
z = np.zeros((ni, nj, nk))

# Accumulate center positions along each axis
for i in range(ni):
    if i == 0:
        x[i, :, :] = 0.5 * di[i, :, :]
    else:
        x[i, :, :] = x[i-1, :, :] + 0.5 * (di[i-1, :, :] + di[i, :, :])

for j in range(nj):
    if j == 0:
        y[:, j, :] = 0.5 * dj[:, j, :]
    else:
        y[:, j, :] = y[:, j-1, :] + 0.5 * (dj[:, j-1, :] + dj[:, j, :])

for k in range(nk):
    if k == 0:
        z[:, :, k] = 0.5 * dk[:, :, k] + depth
    else:
        z[:, :, k] = z[:, :, k-1] + 0.5 * (dk[:, :, k-1] + dk[:, :, k])

#
#
# Read Temperature
#
#

patterns = {
    "TIME":    re.compile(r"TIME\s+=\s+(\d+)"),
    "TEMP":    re.compile(r"^\s*\*?TEMP\s+ALL\s*"),
}

data = { }
prop = []

state = "START"
curr_time_days = 0
with open(temp_fn, 'r') as file:
    ln = 0
    for line in file:
        ln += 1
        line = line.strip()

        if m := patterns["TIME"].search(line) :
            state = "START"
            curr_time_days = float(m.group(1))
            data[curr_time_days] = {}

        line = re.sub(r'\s*\*\*.*', '', line)

        if m := patterns["TEMP"].search(line) :
            state = "READING"
            prop = []
            data[curr_time_days]["TEMP_C"] = prop
            continue

        if state == "READING" :
            ll = _to_list( line )
            prop.extend(ll)

#
#  Temperature for each cell
#
from itertools import product
records = []
for t_days, entry in data.items():
    t_sec = t_days * 24 * 60 * 60
    temp_vals = np.array(entry['TEMP_C']) + 273
    ii = 0
    for k, j, i in product(range(nk), range(nj), range(ni)):
        records.append((t_sec, i, j, k, temp_vals[ii]))
        ii += 1

df_temp = pd.DataFrame(records, columns=['t', 'i', 'j', 'k', 'TEMP'])
df_temp.set_index(['t', 'i', 'j', 'k'], inplace=True)

print(df_temp)



#
#  Create dataframe
#
import pandas as pd
from itertools import product  # All combinations
data = [    (i, j, k, x[i, j, k], z[i, j, k])
            for k, j, i in product(range(nk), range(nj), range(ni)) ]

df = pd.DataFrame(data, columns=['i', 'j', 'k', 'x', 'z'])
df.set_index(['i', 'j', 'k'], inplace=True)

#
# Merge both dataframes
#
df_coords = df.reset_index()
df_temp_reset = df_temp.reset_index()
df_merged = pd.merge(df_temp_reset, df_coords, on=['i', 'j', 'k'])

df = df_merged.set_index(['t', 'i', 'j', 'k'])


import matplotlib.pyplot as plt
from scipy.interpolate import griddata


from matplotlib.ticker import FuncFormatter

# Custom formatter: divide x labels by 10
xscale = 10
def scale_x_ticks(x, _): return f"{x / xscale:.1f}"
def format_time_days(t_days):
    total_hours = int(t_days * 24)
    years, rem_hours = divmod(total_hours, 365 * 24)
    months, rem_hours = divmod(rem_hours, 30 * 24)
    days, hours = divmod(rem_hours, 24)
    return f"{years}y {months}m {days}d {hours}h"

df.to_pickle( "db/temp.pkl" )
df.reset_index()[["t","x","z","TEMP"]].to_csv("db/temperature.csv.gz", index=False, float_format="%.6e", compression="gzip")
print(df)

# #
# # Heat map for each time
# #
# tempmin = df["TEMP"].min()
# tempmax = df["TEMP"].max()
# for t in df.index.get_level_values('t').unique():
#     print(f"Plotting {format_time_days(t)} ... ")
#     df_t = df.loc[t].reset_index()
#     print(df_t)
#     plt.figure()
#     x_vals = df_t['x'].values * xscale
#     z_vals = df_t['z'].values
#     v_vals = df_t["TEMP"].values
#     xi = np.linspace(x_vals.min(), x_vals.max(), 10000)
#     zi = np.linspace(z_vals.min(), z_vals.max(), 300)
#     xi, zi = np.meshgrid(xi, zi)
#     vi = griddata((x_vals, z_vals), v_vals, (xi, zi), method='linear', fill_value=0)
#     plt.imshow(vi, extent=[x_vals.min(), x_vals.max(), z_vals.min(), z_vals.max()],
#                origin='lower', aspect='auto', cmap='coolwarm', interpolation='bilinear',
#                vmin=tempmin, vmax=tempmax)
#     plt.colorbar(label='Value')
#     plt.xlabel('$R$ (m)')
#     plt.ylabel('$TVDSS$ (m)')
#     plt.xlim([0,200*xscale])
#     plt.gca().xaxis.set_major_formatter(FuncFormatter(scale_x_ticks))
#     plt.title(f'Interpolated Heatmap (t={format_time_days(t)})')
#     plt.gca().invert_yaxis()
#     plt.tight_layout()
#     plt.savefig(f'png/temperature_t{t}.png', dpi=150)
#     plt.close()


# # # #
# # # # Scatter
# # # #
# # for t in df.index.get_level_values('t').unique():
# #     df_t = df.loc[t].reset_index()
# #     plt.figure()
# #     x_vals = df_t['x'].values * xscale
# #     z_vals = df_t['z'].values
# #     v_vals = df_t["TEMP"].values
# #     plt.scatter(x_vals, z_vals, c=v_vals, s=1)
# #     plt.xlabel('x')
# #     plt.ylabel('z')
# #     plt.title(f'Scatter (t={format_time_days(t)})')
# #     plt.xlim([0,200*xscale])
# #     plt.gca().xaxis.set_major_formatter(FuncFormatter(scale_x_ticks))
# #     plt.grid(True)
# #     plt.gca().invert_yaxis()
# #     plt.tight_layout()

# #     plt.savefig(f'png/scatter_temperature_t{t}.png', dpi=150)
# #    plt.close()
