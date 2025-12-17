import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

import numpy as np



# Create custom colormap
def create_custom_cmap( vmin, vzero, vmax ):
#     zero = ( vzero - vmin )  / ( vmax - vmin )
    zero = (np.log(vzero) - np.log(vmin)) / (np.log(vmax) - np.log(vmin))

    jet = plt.cm.jet
    colors = [
        (0.0, jet(0.1)),       # -2: blue from jet
        (zero/2, jet(0.3)),
        (zero, 'white'),        # 0: white
#         (zero+(1-zero)/2, jet(0.7)),   ## Enables intermediate yellow
        (1.0, jet(0.92))        # 8: red from jet
    ]
    return mcolors.LinearSegmentedColormap.from_list('custom', colors)

def create_custom_cmap_nozero():
    jet = plt.cm.jet
    colors = [
#         (0.0, jet(0.1)),  
#         (1.0, jet(0.92)) 
#         (0.0, '#6BA3D0'),    # soft blue
        (0.0, '#F5B7B1'),     # soft red
        (1.0, '#FF0000')     # soft red
    ]
    return mcolors.LinearSegmentedColormap.from_list('custom', colors)

def create_custom_cmap_nomin():
    jet = plt.cm.jet
    colors = [
#         (0.0, jet(0.1)),  
#         (1.0, jet(0.92)) 
        (0.0, 'white'),    # soft blue
        (1.0, jet(0.92)) 
#         (1.0, '#E74C3C')     # soft red
    ]
    return mcolors.LinearSegmentedColormap.from_list('custom', colors)

