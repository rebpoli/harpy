#!/usr/bin/env -S python -i

import xarray as xr
import numpy as np

#
#
#
def flip_eigenvector(vec):
    idx = np.argmax(np.abs(vec)) 
    if vec[idx] < 0:
        return -vec
    return vec

#
#
#
def principal_stresses( sigma ) :
    eigvals, eigvecs = np.linalg.eig(sigma)

    # Sort order
    idx = np.argsort(eigvals)
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]    
    
    S1 = eigvals[2]
    S2 = eigvals[1]
    S3 = eigvals[0]
    V1 = flip_eigenvector(eigvecs[:,2])
    V2 = flip_eigenvector(eigvecs[:,1])
    V3 = flip_eigenvector(eigvecs[:,0])

    norm_S1 = np.linalg.norm(V1)
    norm_S2 = np.linalg.norm(V2)
    norm_S3 = np.linalg.norm(V3)

    print( sigma )
    print(f"S1:{S1} => {V1} ({norm_S1})")
    print(f"S2:{S2} => {V2} ({norm_S2})")
    print(f"S3:{S3} => {V3} ({norm_S3})")

    return ( S1, S2, S3, V1, V2, V3 )


#
#
#

ds = xr.open_dataset("points_libmesh.nc")

ds.info()
print("====")
print(ds)

# sigma_ds = ds['sigma']  # shape: (time, x, y, z, 9)

# nt, nx, ny, nz, ncomp = sigma_ds.shape
# assert ncomp == 9, "Expected 9 components for 3x3 tensor"
# # Iterate over all times and points
# for t in range(nt):
#     for x in range(nx):
#         for y in range(ny):
#             for z in range(nz):
#                 # Get the flattened 3x3 tensor
#                 tensor_flat = sigma_ds[t, x, y, z, :].values
#                 # Reshape into 3x3 matrix
#                 sigma = tensor_flat.reshape(3, 3)

#                 print(f"time={t}, x={x}, y={y}, z={z}")
#                 principal_stresses( sigma )

