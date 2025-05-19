#!/usr/bin/env -S python -i

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os, time
from scipy.interpolate import make_interp_spline
from matplotlib.ticker import MultipleLocator
from numpy import exp
from math import log10, log

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

    eps0 = 0.0175 # /s
    sig0 = 9.91e6
    N = 7.55

    sig_bar = sig/sig0

    ## Transient setup
    kappa = 5e-3
    c     = 3e-3
    m     = 1 
    eps_tr_star = kappa * exp( c * T ) * pow( sig_bar, m )

    alpha_w = 2.4
    beta_w =  0
    alpha_r = 0
    beta_r = 0
    print(f"eps^* = {eps_tr_star}")

    eps_ss_rate = eps0 * exp( -Q/R/T ) * sig_bar**N

    eps_tot = np.zeros_like(t_s)
    eps_tr = np.zeros_like(t_s)
    F_i = np.zeros_like(t_s)
    Zeta_i = np.zeros_like(t_s)

    for i in range( len( t_s ) ) :
        ess_n = 0 ; dt = 0; etr_n=0
        if i > 0 : 
            ess_n = eps_tot[i-1]
            etr_n = eps_tr[i-1]
            dt = t_s[i] - t_s[i-1]

        t = t_s[i]

        zeta = 1 - etr_n/eps_tr_star

        alpha = alpha_w
        beta = beta_w 
        if zeta < 0 : 
            alpha=-alpha_r
            beta=-beta_r 

        test = alpha + beta * log10(sig_bar)
#         print(f"test:{test:.2e} alpha:{alpha:.2e} beta:{beta:.2e} sig_bar:{sig_bar:.2e} zeta:{zeta:.2e}")
        F = exp(alpha*zeta*zeta) * pow( sig_bar, beta*zeta*zeta )
        if test < 0 : F = 1

        eps_tr_rate = ( F - 1 ) * eps_ss_rate
        eps_tr[i] = etr_n + eps_tr_rate * dt

        print(f"F:{F:.2e} eps_tr_rate:{eps_tr_rate:.2e} eps_tr:{eps_tr[i]:.2e} etr_n:{etr_n:.2e}")

        eps_tot[i] = ess_n + F * eps_ss_rate * dt

        ## Debugging stuff
        F_i[i] = F
        Zeta_i[i] = zeta


    # Do the plotting
    fig, [ax1, ax3] = plt.subplots(2,1,figsize=(10, 6))
    ax2 = ax1.twinx()

    lines = []
    labels = []


    l1 = ax1.scatter(df.time_h, df.epsxx, c='b', s=10, label=r'$\varepsilon_{xx}$') 
    l2 = ax2.scatter(df.time_h, df.eps_rate, c='r', s=10, label=r'$\dot{\varepsilon}$') 
    l3, = ax2.plot(df.time_h, df.eps_rate_smooth, 'r--') 

    l4, = ax1.plot( t_h, eps_tot , 'k--', lw=1.5, alpha=0.6, label=r"Model ($\varepsilon_{tr}$)" ) 
    l5, = ax1.plot(t_h, eps_tr, 'g--', lw=1.5, alpha=0.6, label=r"Model ($\varepsilon_{tr}$)") 

    lines = [ l1, l2, l4, l5 ]
    labels = [line.get_label() for line in lines]
    ax1.legend(lines, labels, loc='upper left', bbox_to_anchor=(0.04,0.96), fontsize=10 )

    set_grid( ax1, ax2 )

    # Aux variables plot
    l1, = ax3.plot(t_h, F_i, label="$F$" )
    ax3.set_ylabel("$F$")

    ax4 = ax3.twinx()
    l2, = ax4.plot(t_h, Zeta_i, label=r"$\zeta$" , c='red' )
    ax4.set_ylabel(r"$\zeta$")

    lines = [ l1, l2 ]
    labels = [line.get_label() for line in lines]
    ax3.legend(lines, labels, loc='upper left', bbox_to_anchor=(0.3,0.96), fontsize=10 )

    ax3.set_xlim( 0 , ax1.get_xlim()[1])
    ax4.grid(False, which='both')

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
