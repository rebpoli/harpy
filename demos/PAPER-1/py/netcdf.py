import xarray as xr

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


def read_netcdf( fn ) :
    ds = xr.open_dataset(fn)
    ds = assign_component_names(ds)
    return ds
