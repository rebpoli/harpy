
import xarray as xr
import numpy as np
import pandas as pd
import os

## Assign component labels for vec3_comp and ten9_comp based on attrs['components']
def assign_component_names(ds):
    updated = ds
    for var in ds.data_vars:
        if "components" in ds[var].attrs:
            comps = ds[var].attrs["components"].split()

            if "vec3_comp" in ds[var].dims and len(comps) == len(ds["vec3_comp"]):
                updated = updated.assign_coords(vec3_comp=comps)

            if "ten9_comp" in ds[var].dims and len(comps) == len(ds["ten9_comp"]):
                updated = updated.assign_coords(ten9_comp=comps)
    return updated

#
#
#
def read_netcdf( fn ) :
    ds = xr.open_dataset(fn)
    ds = assign_component_names(ds)
    ds['time_in_days'] = ds['time'] / 86400
    return ds

#
#
#
def filter_time( ds, target_interval, max_time ) :
    time_values = ds.time.values
    selected_indices = [0]  
    current_time = time_values[0]

    for i in range(1, len(time_values)):
        # Calculate time difference
        time_diff = time_values[i] - current_time

        # Check if at least 1 hour has passed since last selected point
        if time_diff >= target_interval:
            selected_indices.append(i)
            current_time = time_values[i]

            if (current_time - time_values[0]) >= max_time : break

    return ds.isel(time=selected_indices)

#
#
#
def select_nearest_point(ds, pt, quiet=0):
    if not quiet : print(f"Selecting nearest point to {pt} ...")
    distances = []
    for i in range(len(ds.point_idx)):
        dist_sq = 0
        for j in range(3) :
            c = ds['Coord'].isel(point_idx=i, vec3_comp=j).values
            dist_sq += (c - pt[j]) ** 2

        dist = np.sqrt(dist_sq)
        distances.append( dist );

    nearest_idx = np.argmin(distances)
    spt = ds['Coord'].isel(point_idx=nearest_idx).values.tolist()

    print(f"Selected point: {spt}")
    return ds.isel(point_idx=nearest_idx)


#
#
#
def cache_outdated(source_file, cache_file):
    if not os.path.exists(cache_file): return True
    if not os.path.exists(source_file): return False

    source_time = os.path.getmtime(source_file)
    cache_time = os.path.getmtime(cache_file)

    return source_time > cache_time

#
#
#
def cached( builder_func, cache_file, source_file, force=0 ):
    basedir = "run/cache/"
    cache_file = basedir + cache_file

    needs_rebuild = cache_outdated( source_file, cache_file )
    needs_rebuild |= force

    if needs_rebuild:
        print(f"Building dataset and saving to {cache_file}...")
        ds = builder_func()

        os.makedirs(basedir, exist_ok=True)
        ds.to_netcdf(cache_file)
        return ds
    else:
        print(f"Loading from cache: {cache_file}...")
        return xr.open_dataset(cache_file)

#
#
#
def extract_dataframe(ds, vnames):
    if isinstance(vnames, str): vnames = [vnames]
    import re

    # Name maps
    simple_var_mapping = {
        'P': 'Pressure',
        'T': 'Temperature',
        'DELTAP': 'Delta_P',
        'DELTAT': 'Delta_T',
    }

    vec3_vars = ['U', 'S1', 'S2', 'S3']

    ten9_var_mapping = {
        'SIGTOT': 'Total Stress',
        'SIGEFF': 'Effective Stress',
        'SIGE': 'Elastic Stress',
        'SIGV': 'Viscous Stress',
        'STRAINVP': 'VP Strain',
        'STRAINRATEVP': 'VP Strain Rate',
    }

    # Regular expressions for pattern matching
    vec3_pattern = re.compile(r'^(.+)([XYZ])$', re.IGNORECASE)
    ten9_pattern = re.compile(r'^(.+)(XX|XY|XZ|YX|YY|YZ|ZX|ZY|ZZ)$', re.IGNORECASE)

    data_dict = { 'time_in_days': ds['time_in_days'].values }

    for vname in vnames:
        # Check simple variable mapping
        if vname in simple_var_mapping:
            ds_var = simple_var_mapping[vname]
            data_dict[vname] = ds[ds_var].values

        elif (match := ten9_pattern.match(vname)):
            base_name = match.group(1)
            comp = match.group(2).lower()

            if base_name in ten9_var_mapping:
                ds_var = ten9_var_mapping[base_name]
                data_dict[vname] = ds[ds_var].sel(ten9_comp=comp).values

        elif (match := vec3_pattern.match(vname)):
            base_name = match.group(1)
            comp = match.group(2).lower()

            if base_name in vec3_vars:
                data_dict[vname] = ds[base_name].sel(vec3_comp=comp).values

        elif vname in ds:
            data_dict[vname] = ds[vname].values

        else:
            print(f"Warning: Variable '{vname}' not found or not recognized")

    # Create DataFrame with time as index
    df = pd.DataFrame( data_dict )

    return df
