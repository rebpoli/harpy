#!/usr/bin/env -S python -i
import sys
sys.path.append('../py')
sys.stdout.reconfigure(line_buffering=True)


from netcdf import read_netcdf, filter_time, select_nearest_point, cached
import numpy as np

import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')  ## Noninteractive

#
# Dataset transformation routines
#
def _reduce() :
    global ds
    print("Reducing dataset ...")
    target_interval = 1 * 60 * 60   # 1h => 3600s
    max_time = 100 * 24 * 60 * 60   # 100 days => seconds
    return filter_time( ds , target_interval , max_time )
#
def _select_point() :
    global ds_reduced
    print("Selecting point ...")
    ds_point = select_nearest_point( ds_reduced, [5,5,5] )
    return ds_point

## Dataset manipulation
src_cdf = "run/cdf/lin_0.cd" 
ds = read_netcdf(src_cdf)
ds_reduced = cached( _reduce, "ds_reduced.cd", src_cdf )
ds_point = cached( _select_point , "ds_point.cd", src_cdf, 0 )

ds_point.load()
uy = ds_point['U'].sel(vec3_comp='y')

plt.figure(figsize=(10, 6))
plt.plot(ds_point.time_in_days.values, uy.values)
plt.xlabel('Time')
plt.ylabel('UY')
plt.grid(True)
plt.tight_layout()

plt.savefig('uy.png', dpi=300, bbox_inches='tight')
plt.close()

