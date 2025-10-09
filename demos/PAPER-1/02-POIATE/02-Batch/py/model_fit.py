#!/usr/bin/env -S python -i

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os, time
from scipy.interpolate import make_interp_spline
from matplotlib.ticker import MultipleLocator
from numpy import exp
from math import log10, log

#
#
#
def build_model_df( t_s , T_C, sig_MPa ) :

    sig = sig_MPa * 1e6
    T = T_C + 273
    fud = 1.08  ##  A Fudge factor (see arrhenius)

    Q = 51600
    R = 8.314

    ## Steady State
    sig0ss = [ 10e6 ,   50e6  ]
    N    = [    8.6 ,   .7    ]

    ## Transient setup
    sig0tr     = 3e9
    c          = 3e-3
    m          = 1
    alpha_w    = 2


    # Steady State
    eps_ss_rate = 0
    for i in range(2) : 
        eps_ss_rate += exp( -(Q/R/T)**fud  ) * ( sig/sig0ss[i] )**N[i]

    eps_tr_star = exp( c * T ) * pow( sig/sig0tr, m )

    eps_tot = np.zeros_like(t_s)
    eps_tr = np.zeros_like(t_s)
    F_i = np.zeros_like(t_s)
    Zeta_i = np.zeros_like(t_s)

    # Integrate in time
    for ti in range( len( t_s ) ) :
        ess_n = 0 ; dt = 0; etr_n=0
        if ti > 0 : 
            ess_n = eps_tot[ti-1]
            etr_n = eps_tr[ti-1]
            dt = t_s[ti] - t_s[ti-1]

        t = t_s[ti]

        zeta = 1 - etr_n/eps_tr_star

        alpha = alpha_w
        if zeta < 0 :  alpha = 0 

        F = exp(alpha*zeta*zeta)

        eps_tr_rate = ( F - 1 ) * eps_ss_rate
        eps_tr[ti] =  etr_n + eps_tr_rate * dt
        eps_tot[ti] =  ess_n + F * eps_ss_rate * dt

        ## Debugging stuff
        F_i[ti] = F
        Zeta_i[ti] = zeta

    ## Feed the dataframe
    df = pd.DataFrame( { 
                         't_s'           : t_s,
                         't_h'           : t_s/60/60,
                         'T'             : T,
                         'T_C'           : T_C,
                         'sig'           : sig,
                         'sig_MPa'       : sig/1e6,
                         'eps_tot'       : eps_tot,
                         'eps_ss_rate'   : eps_ss_rate,
                         'eps_ss_rate_h' : eps_ss_rate*60*60,
                         'eps_tr'        : eps_tr,
                         'F'             : F_i ,
                         'Zeta'          : Zeta_i
                        } )
    return df

import itertools
colors = ['red', 'blue', 'green', 'orange', 'purple', 'yellow', 'pink', 'brown', 'gray', 'cyan', 'magenta', 'lime', 'indigo', 'violet', 'turquoise']

#
#
#
def run() :
    fn = "py/paper.mplstyle"
    if os.path.isfile(fn) :
        plt.style.use('default')   ## reset!
        plt.style.use(fn)

    # Do the plotting
    fig, [ [ ax1, ax2 ] , [ ax3, ax4 ] ] = plt.subplots(2,2,figsize=(10, 6))
    time_plots( ax1, ax2 )
    sige_plot( ax3 )

    plt.show()
#
#
#
def sige_plot( ax ) :
    t_h = np.array( [ 24 * 3600 ] ) # 1d


    dfs = []
    for T_C in [ 130, 86, 43 ] :
        for SIGE in 10**np.arange(0.1, 1.7, 0.1) :
            dfs.append( build_model_df( t_h, T_C, SIGE ) )
    df = pd.concat(dfs)

    color_cycle = itertools.cycle(colors)
    for T_C, gdf in df.groupby( 'T_C' ) :
        c = next(color_cycle)
        ax.plot( gdf.sig_MPa, gdf.eps_ss_rate, c=c, label=f"{T_C} C" )
#         ax.plot( gdf.sig_MPa, gdf.eps_ss_rate1_h, ':',c=c )

    fac = 1.08
    Q = 51600
    R = 8.314

    color_cycle = itertools.cycle(colors)
    df = pd.read_csv("py/csv/fig_5_17_T43.csv", sep="\t")
    l1 = ax.scatter(df.sig_MPa, df.eps_ss_rate_h/60/60 , c=next(color_cycle), s=10, label=r'$\varepsilon_{xx}$ $T=43$ °C') 
    df = pd.read_csv("py/csv/fig_5_17_T86.csv", sep="\t")
    l1 = ax.scatter(df.sig_MPa, df.eps_ss_rate_h/60/60 , c=next(color_cycle), s=10, label=r'$\varepsilon_{xx}$ $T=86$ °C') 
    df = pd.read_csv("py/csv/fig_5_17_T130.csv", sep="\t")
    l1 = ax.scatter(df.sig_MPa, df.eps_ss_rate_h/60/60 , c=next(color_cycle), s=10, label=r'$\varepsilon_{xx}$ $T=130$ °C') 

    ax.legend()
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel(r"$\sigma_e$ (MPa)")
    ax.set_ylabel(r"$\dot\varepsilon$ (1/s)")

#
#
#
def time_plots( ax1, ax2 ) :
    # Solve the model
    t_h = np.arange( 0, 2200, 10) 
    t_s = t_h * 60. * 60.

    # Compute the model
    dfs_16mpa = []
    for T_C in [ 16, 36, 56 , 86 ] :
        dfs_16mpa.append ( build_model_df( t_s, T_C, 16 ) )
    df_16mpa = pd.concat( dfs_16mpa )

    dfs_26mpa = []
    for T_C in [ 16, 36, 56, 86 ] :
        dfs_26mpa.append ( build_model_df( t_s, T_C, 26e6 ) )
    df_26mpa = pd.concat( dfs_26mpa )


    ax1_lines, ax1_labels = [], []
    ax2_lines, ax2_labels = [], []

    # 16 MPa
    color_cycle = itertools.cycle(colors)
    for T_C, gdf in df_16mpa.groupby( 'T_C' ) :
        t_lab = f"$T={T_C}$ °C"
        c = next(color_cycle)
        l1, = ax1.plot( gdf.t_h, gdf.eps_tot , '--', lw=1.5, alpha=0.6, label=r"Model ($\varepsilon_{tot}$) "+t_lab, c=c ) 
        ax1_lines.append(l1)
        ax1_labels.append(l1.get_label())
    ax1.set_title(r"$\sigma_e= 16$ MPa")
    ax1.legend(ax1_lines, ax1_labels, loc='upper left', bbox_to_anchor=(0.04,0.96), fontsize=10 )

    # 26 MPa
    color_cycle = itertools.cycle(colors)
    for T_C, gdf in df_26mpa.groupby( 'T_C' ) :
        c = next(color_cycle)
        l1, = ax2.plot( gdf.t_h, gdf.eps_tot , '--', lw=1.5, alpha=0.6, label=r"Model ($\varepsilon_{tot}$) "+t_lab, c=c ) 
        ax2_lines.append(l1)
        ax2_labels.append(l1.get_label())
    ax2.set_title(r"$\sigma_e= 26$ MPa")
    ax2.legend(ax2_lines, ax2_labels, loc='upper left', bbox_to_anchor=(0.04,0.96), fontsize=10 )

    for ax in [ ax1, ax2 ] :
        ax.set_ylabel(r"strain  $\varepsilon$ (-)")
        ax.set_xlabel(r"time (h)")

#     #
#     #  Data from Poiate, 2012
#     #
    df = pd.read_csv("py/csv/fig_5_15.csv", sep="\t")
#     df["eps_rate"] = df.epsxx.diff() / df.time_h.diff()
#     df['eps_rate'] = df['eps_rate'].replace([np.inf, -np.inf], None)
#     df['eps_rate'] = df['eps_rate'].where(pd.notna(df['eps_rate']), None)
#     df['eps_rate_smooth'] = df['eps_rate'].rolling(window=2, center=True).mean()
    l1 = ax1.scatter(df.time_h, df.epsxx, c='b', s=10, label=r'$\varepsilon_{xx}$') 
#     ax1_lines.extend( [ l1, l2 ] )
#     ax1_labels.extend( [line.get_label() for line in [l1,l2]] )

    #



#
#
#
def set_grid( ax1, ax1b ) :
    ax1.set_ylim( 0 , ax1.get_ylim()[1])
    ax1.set_xlim( 0 , ax1.get_xlim()[1])

    ax1.grid(True)
    ax1b.grid(False, which='both')

    from matplotlib.ticker import ScalarFormatter
    formatter = ScalarFormatter(useMathText=True)
    formatter.set_powerlimits((-4, -4))  # Force scientific notation with 10^-4
    ax1b.yaxis.set_major_formatter(formatter)

    y1_range = ax1.get_ylim()[1] - ax1.get_ylim()[0]
    y2_range = ax1b.get_ylim()[1] - ax1b.get_ylim()[0]
    ratio = y2_range / y1_range    

    y2_min, y2_max = ax1b.get_ylim()
    y1_min, y1_max = ax1.get_ylim()

    ax1.set_ylabel(r"strain  $\varepsilon$ (mm/mm)")
    ax1b.set_ylabel(r"strain rate  $\dot\varepsilon$ (/h)")
    ax1.set_xlabel(r"time (h)")

    y0 = 5e-5

    ax1_ticks = ax1.get_yticks()
    ax1b.set_yticks(ax1_ticks * ratio + y0)
    ax1b.set_ylim(y1_min * ratio + y0 , y1_max * ratio + y0 )


####
run()
