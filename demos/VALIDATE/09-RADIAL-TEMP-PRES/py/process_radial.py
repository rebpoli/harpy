#!/usr/bin/env -S python -i

import os, time
import re
import numpy as np
import pandas as pd
import gzip

PROP = "Pressure"
# PROP = "Temperature"

# Parse the dat to capture the radial grid
dat_fn = "stars/test3-refine-dates-expon.dat"
temp_fn = f"stars/test3-refine-dates-expon {PROP}.txt.gz"

# Precompile regex patterns for performance
patterns = {
    "GRID":    re.compile(r"^\s*\*?GRID\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+([0-9.]+)$"),
    "DI":      re.compile(r"^\s*\*?DI\s+(\S+)\s+(.+)?$"),
    "DJ":      re.compile(r"^\s*\*?DJ\s+(\S+)\s+(.+)?$"),
    "DK":      re.compile(r"^\s*\*?DK\s+(\S+)\s+(.+)?$"),
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
with open(dat_fn, 'r', encoding="latin1") as file:
    print(f"Openin ... {dat_fn}")
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
            print(f"ni:{ni} nj:{nj} nk:{nk} rwstr:{rwstr} rw:{rw}")
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

print(di)

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
    "TIME":    re.compile(r"TIME\s+=\s+(\d+)\s+(\S+)(?:\s+(\S+))?"),
    "TEMP":    re.compile(r"^\s*\*?TEMP\s+ALL\s*"),
    "PRES":    re.compile(r"^\s*\*?PRES\s+ALL\s*"),
}

##
#
#
def fraction_of_day(curr_time) :
    from datetime import datetime, timedelta
    t = datetime.strptime(curr_time, "%H:%M:%S")
    td = timedelta(hours=t.hour, minutes=t.minute, seconds=t.second)
    return round( td.total_seconds() / 86400 , 6 )


data = { }
prop = []

state = "START"
curr_time_days = 0
with gzip.open(temp_fn, 'rt') as file:
    ln = 0
    for line in file:
        ln += 1
        line = line.strip()

        if m := patterns["TIME"].search(line) :
            state = "START"
            curr_time_days = float(m.group(1))
            curr_date = m.group(2)
            curr_time = m.group(3)
            if not curr_time : curr_time = "00:00:00"

            curr_time_days += fraction_of_day( curr_time )
            data[curr_time_days] = {}
            print(f"Reading time {curr_time_days} ({curr_date} {curr_time})")

        line = re.sub(r'\s*\*\*.*', '', line)

        if m := patterns["TEMP"].search(line) :
            state = "READING"
            prop = []
            data[curr_time_days]["TEMP_C"] = prop
            continue

        if m := patterns["PRES"].search(line) :
            state = "READING"
            prop = []
            data[curr_time_days]["PRES_KGF"] = prop
            continue

        if state == "READING" :
            ll = _to_list( line )
            prop.extend(ll)

#
#  Each cell Temp or Pres
#

print("Building dataframe ...")
from itertools import product

if PROP == "Temperature" :
    records = []
    for t_days, entry in data.items():
        t_sec = t_days * 24 * 60 * 60
        temp_vals = np.array(entry['TEMP_C']) + 273
        ii = 0
        for k, j, i in product(range(nk), range(nj), range(ni)):
            if len(temp_vals) :
                records.append((t_sec, i, j, k, temp_vals[ii]))
            ii += 1

    df_temp = pd.DataFrame(records, columns=['t', 'i', 'j', 'k', 'TEMP'])
    df_temp.set_index(['t', 'i', 'j', 'k'], inplace=True)

    print(df_temp.index.get_level_values("t"))

if PROP == "Pressure" :
    records = []
    for t_days, entry in data.items():
        t_sec = t_days * 24 * 60 * 60
        pres_vals = np.array(entry['PRES_KGF']) * 98066.5
        ii = 0
        for k, j, i in product(range(nk), range(nj), range(ni)):
            if len(pres_vals) :
                records.append((t_sec, i, j, k, pres_vals[ii]))
            ii += 1

    df_temp = pd.DataFrame(records, columns=['t', 'i', 'j', 'k', 'PRES'])
    df_temp.set_index(['t', 'i', 'j', 'k'], inplace=True)

    print(df_temp.index.get_level_values("t"))


#
#  Create dataframe
#
print("Organizing dataframe ...")
import pandas as pd
from itertools import product  # All combinations
data = [    (i, j, k, x[i, j, k], z[i, j, k])
            for k, j, i in product(range(nk), range(nj), range(ni)) ]

df = pd.DataFrame(data, columns=['i', 'j', 'k', 'x', 'z'])
df.set_index(['i', 'j', 'k'], inplace=True)

#
# Merge both dataframes
#
print("Merging dataframe ...")
df_coords = df.reset_index()
df_temp_reset = df_temp.reset_index()
df_merged = pd.merge(df_temp_reset, df_coords, on=['i', 'j', 'k'])

df = df_merged.set_index(['t', 'i', 'j', 'k'])


import matplotlib.pyplot as plt
from scipy.interpolate import griddata


from matplotlib.ticker import FuncFormatter

print(df.index.get_level_values("t"))

print("Exporting pickle ...")
# df.to_pickle( "model/temp.pkl" )
# print("Exporting gzip ...")
if PROP == "Temperature" : field = "TEMP"
if PROP == "Pressure" : field = "PRES"

df.reset_index()[["t","x","z",field]].to_csv(f"model/{PROP.lower()}.csv.gz", index=False, float_format="%.6e", compression="gzip")

