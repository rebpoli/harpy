import xarray as xr
import numpy as np

# Return indices of points that are at least min_dist apart.
def spatial_subsample(Y, Z, min_dist):
    keep_idx = []
    for i, (y, z) in enumerate(zip(Y, Z)):
        if all((y - Y[j])**2 + (z - Z[j])**2 >= min_dist**2 for j in keep_idx):
            keep_idx.append(i)
    return np.array(keep_idx)

