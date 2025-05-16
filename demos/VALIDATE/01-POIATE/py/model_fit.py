#!/usr/bin/env -S python -i

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os, time
from scipy.interpolate import make_interp_spline
from matplotlib.ticker import MultipleLocator

def run() :
    fn = "py/paper.mplstyle"
    if os.path.isfile(fn) :
        plt.style.use('default')   ## reset!
        plt.style.use(fn)


    df = pd.read_csv("py/csv/fig15_5.csv", sep="\t")

    df["eps_rate"] = df.epsxx.diff() / df.time_h.diff()
    df['eps_rate'] = df['eps_rate'].replace([np.inf, -np.inf], None)
    df['eps_rate'] = df['eps_rate'].where(pd.notna(df['eps_rate']), None)

    df['eps_rate_smooth'] = df['eps_rate'].rolling(window=2, center=True).mean()

    # Solve the model
    t_h = np.arange( 0, 2200, 10) 
    t_s = t_h * 60. * 60.

    Q = 51600
    R = 8.314
    T = 86 + 273
    sig = 16e6

    eps0 = 0.017 # /s
    sig0 = 9.91e6
    N = 7.55

    sig_bar = sig/sig0

    eps_ss_rate = eps0 * np.exp( -Q/R/T ) * sig_bar**N

    eps_tot = np.zeros_like(t_s)
    eps_tr = np.zeros_like(t_s)

    for i in range( len( t_s ) ) :
        en = 0 ; dt = 0
        if i > 0 : 
            en = eps_tot[i-1]
            dt = t_s[i] - t_s[i-1]

        t = t_s[i]
        eps_tot[i] = en + eps_ss_rate * dt


    # Do the plotting
    fig, ax1 = plt.subplots(figsize=(10, 6))
    ax1.scatter(df.time_h, df.epsxx, c='b', s=10, label=r'$\varepsilon_{xx}$')

#     ax1.plot( t_h, eps_tot , 'k--' )

    ax2 = ax1.twinx()
    ax2.scatter(df.time_h, df.eps_rate, c='r', s=10, label=r'$\dot{\varepsilon}$')
    ax2.plot(df.time_h, df.eps_rate_smooth, 'r--')

    set_grid( ax1, ax2 )

    fig.legend( loc='upper left', bbox_to_anchor=(0.07,0.96), fontsize=14 )

    plt.show()


#
#
#
def set_grid( ax1, ax2 ) :
    ax1.set_ylim( 0 , ax1.get_ylim()[1])
    ax1.set_xlim( 0 , ax1.get_xlim()[1])

    ax1.grid(True)
    ax2.grid(False, which='both')

    from matplotlib.ticker import ScalarFormatter
    formatter = ScalarFormatter(useMathText=True)
    formatter.set_powerlimits((-4, -4))  # Force scientific notation with 10^-4
    ax2.yaxis.set_major_formatter(formatter)

    y1_range = ax1.get_ylim()[1] - ax1.get_ylim()[0]
    y2_range = ax2.get_ylim()[1] - ax2.get_ylim()[0]
    ratio = y2_range / y1_range    

    y2_min, y2_max = ax2.get_ylim()
    y1_min, y1_max = ax1.get_ylim()

    ax1.set_ylabel(r"Strain  $\varepsilon$ (mm/mm)")
    ax2.set_ylabel(r"Strain Rate  $\dot\varepsilon$ (/h)")
    ax1.set_xlabel(r"Time (h)")

    y0 = 5e-5

    ax1_ticks = ax1.get_yticks()
    ax2.set_yticks(ax1_ticks * ratio + y0)
    ax2.set_ylim(y1_min * ratio + y0 , y1_max * ratio + y0 )


####
run()
