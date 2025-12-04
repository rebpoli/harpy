import numpy as np
import pandas as pd


#
#
#
def add_noflow_ghosts(df, z0_list):
    ghost_rows = []
    df['is_ghost'] = False

    print(f"Processing {len(z0_list)} interface(s)...")

    for idx, z_0 in enumerate(z0_list):
        print(f"\n[{idx+1}/{len(z0_list)}] Processing interface at z_0 = {z_0}")
        eps = 0.1

        def process_column(group):
            # Sort by z
            group = group.sort_values('z')
            z_vals = group['z'].values
            p_vals = group['PRES'].values

            ghosts = []

            # Positive side: z > z_0
            positive_mask = z_vals > z_0
            if np.sum(positive_mask) >= 2:
                z_positive = z_vals[positive_mask]
                p_positive = p_vals[positive_mask]

                # Find two closest points to z_0 + eps
                distances = np.abs(z_positive - (z_0 + eps))
                closest_indices = np.argsort(distances)[:2]

                z1, z2 = z_positive[closest_indices]
                p1, p2 = p_positive[closest_indices]

                if z2 != z1:
                    a = (p2 - p1) / (z2 - z1)
                    b = p1 - a * z1
                    p_ghost = a * (z_0 + eps) + b

                    ghosts.append({
                        't': group['t'].iloc[0],
                        'x': group['x'].iloc[0],
                        'z': z_0 + eps,
                        'PRES': p_ghost,
                        'is_ghost': True
                    })

            # Negative side: z < z_0
            negative_mask = z_vals < z_0
            if np.sum(negative_mask) >= 2:
                z_negative = z_vals[negative_mask]
                p_negative = p_vals[negative_mask]

                distances = np.abs(z_negative - (z_0 - eps))
                closest_indices = np.argsort(distances)[:2]

                z1, z2 = z_negative[closest_indices]
                p1, p2 = p_negative[closest_indices]

                if z2 != z1:
                    a = (p2 - p1) / (z2 - z1)
                    b = p1 - a * z1
                    p_ghost = a * (z_0 - eps) + b

                    ghosts.append({
                        't': group['t'].iloc[0],
                        'x': group['x'].iloc[0],
                        'z': z_0 - eps,
                        'PRES': p_ghost,
                        'is_ghost': True
                    })

            return ghosts

        # Count total groups for progress tracking
        groups = list(df.groupby(['t', 'x']))
        total_groups = len(groups)
        print(f"  Total (t, x) columns to process: {total_groups}")

        ghosts_added = 0
        # Process all (t, x) columns at once using groupby
        for i, (_, group) in enumerate(groups):
            if (i + 1) % 100 == 0 or (i + 1) == total_groups:
                print(f"  Progress: {i+1}/{total_groups} columns processed ({100*(i+1)/total_groups:.1f}%)", end='\r')

            ghosts = process_column(group)
            ghosts_added += len(ghosts)
            ghost_rows.extend(ghosts)

        print(f"\n  Added {ghosts_added} ghost nodes for this interface")

    # Convert ghost rows to dataframe and concatenate with original
    if ghost_rows:
        ghost_df = pd.DataFrame(ghost_rows)
        df = pd.concat([df, ghost_df], ignore_index=True)
        print(f"\n✓ Total ghost nodes added: {len(ghost_rows)}")
        print(f"✓ Final dataframe size: {len(df)} rows")
    else:
        print("\n⚠ No ghost nodes were added")

    return df


# #
# #
# #
# def add_noflow_ghosts(df, z0_list):
#     ghost_rows = []

#     T = df['t'].values
#     X = df['x'].values
#     Z = df['z'].values
#     P = df['PRES'].values
#     df['is_ghost'] = False

#     for z_0 in z0_list:
#         tolerance = 13
#         mask = np.abs(Z - z_0) < tolerance
#         print(df[mask])

#         if np.any(mask):
#             t_near = T[mask]
#             x_near = X[mask]
#             z_near = Z[mask]
#             p_near = P[mask]

#             eps = 0.1
#             for i in range(len(x_near)):
#                 z_ghost = z_0 - eps if z_near[i] < z_0 else z_0 + eps
#                 ghost_row = { 't':t_near[i], 'x': x_near[i], 'z': z_ghost, 'PRES': p_near[i], 'is_ghost': True }
#                 ghost_rows.append(ghost_row)

#     # Convert ghost rows to dataframe and concatenate with original
#     if ghost_rows:
#         ghost_df = pd.DataFrame(ghost_rows)
#         df = pd.concat([df, ghost_df], ignore_index=True)

#     return df

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
