import numpy as np
import pandas as pd

#
#
#
def add_noflow_ghosts(df, z0_list):
    ghost_rows = []

    T = df['t'].values
    X = df['x'].values
    Z = df['z'].values
    P = df['PRES'].values
    df['is_ghost'] = False

    for z_0 in z0_list:
        tolerance = 13
        mask = np.abs(Z - z_0) < tolerance
        print(df[mask])

        if np.any(mask):
            t_near = T[mask]
            x_near = X[mask]
            z_near = Z[mask]
            p_near = P[mask]

            eps = 0.1
            for i in range(len(x_near)):
                z_ghost = z_0 - eps if z_near[i] < z_0 else z_0 + eps
                ghost_row = { 't':t_near[i], 'x': x_near[i], 'z': z_ghost, 'PRES': p_near[i], 'is_ghost': True }
                ghost_rows.append(ghost_row)

    # Convert ghost rows to dataframe and concatenate with original
    if ghost_rows:
        ghost_df = pd.DataFrame(ghost_rows)
        df = pd.concat([df, ghost_df], ignore_index=True)

    return df

#
#
def get_noflow_ghosts( df, z0_list ) :

    # Add ghost points (same code as before)
    x_bg_aug = df['x'].values.copy()
    z_bg_aug = df['z'].values.copy()
    dt_bg_aug = df['PRES'].values.copy()

    xg = []
    zg = []
    dg = []

    for z_0 in z0_list:
        tolerance = 15
        mask = np.abs(z_bg_aug - z_0) < tolerance

        if np.any(mask):
            x_near = x_bg_aug[mask]
            z_near = z_bg_aug[mask]
            dt_near = dt_bg_aug[mask]
            eps = 1

            for i in range(len(x_near)):
                z_ghost = z_0 - eps if z_near[i] < z_0 else z_0 + eps
                xg.append(x_near[i])
                zg.append(z_ghost)
                dg.append(dt_near[i])

    return xg, zg, dg
