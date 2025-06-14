#!/usr/bin/env -S python -i

import pandas as pd
import matplotlib.pyplot as plt
import os, time
from env import wlog, ilog, dlog, elog, flog
import numpy as np
from glob import glob
from os.path import isdir
import re
import subprocess

def run_main() :
#     fn = "py/paper.mplstyle"
#     if os.path.isfile(fn) :
#         plt.style.use('default')   ## reset!
#         plt.style.use(fn)

    rdirs = [ d for d in glob('batch/run.*/') if isdir(d) ]

    ALL_REG = []
    ALL_SCALARS = []

    for rdir in rdirs :
        ilog(f"Processing {rdir} ...")

        reg = process_slurm_out( rdir )
        ALL_REG.append(reg)

        df = process_scalars(rdir)
        process_model(rdir, df)
        ALL_SCALARS.append(df)


    # Runtime
    df1 = pd.DataFrame( ALL_REG )
    print(df1)

    # eps rate etc
    df = pd.concat( ALL_SCALARS )
    df_ = df.drop_duplicates(["sig", "T", "eps_yy_rate" ])
    print(df_)

    fig, [ax1, ax2] = plt.subplots(1, 2, figsize=(15,6))
    
    SEL_TEMP = [] # np.array([ 130, 86, 46 ]) + 273

    import itertools
    colors = ['red', 'blue', 'green', 'orange', 'purple', 'yellow', 'pink', 'brown', 'gray', 'cyan', 'magenta', 'lime', 'indigo', 'violet', 'turquoise']
    color_cycle = itertools.cycle(colors)

    groups = list(df_.groupby("T", sort=True))
    for k, g in reversed(groups):
        print(k) 
        if len(SEL_TEMP) and not k in SEL_TEMP : continue
        ax1.scatter( -g.sig/1e6, -g.eps_yy_rate, marker='x', lw=1, label=f"T={k-273:.1f} °C", c=next(color_cycle))

    color_cycle = itertools.cycle(colors)
    all_temp = df["T"].drop_duplicates()
    for temp in reversed(list(all_temp)) :
        if len(SEL_TEMP) and not temp in SEL_TEMP : continue
        all_rate = []
        all_sig = np.logspace(np.log10(1e6), np.log10(50e6), num=15)
        for sig in all_sig :
            all_rate.append( analytical_model( temp-273, sig/1e6 ) )
        ax1.plot( all_sig/1e6, all_rate, '-', c=next(color_cycle) ) #, label=f"Analytical $T={temp-273}$ °C", c=next(color_cycle) )

    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.legend()
    ax1.set_ylabel(r"$\dot{\varepsilon}_{yy}$ (1/s)")
    ax1.set_xlabel(r"$\sigma_{yy}$ (MPa)")
    ax1.set_title("Creep strain rate (1/s)")

    ## Plot a colormap
    x = -df["sig"].values / 1e6  # MPa
    y = df["T"].values - 273 # ºC
    z = -df['eps_yy_rate'].values

    # cria a malha
    grid_x, grid_y = np.mgrid[min(x):max(x):100j, min(y):max(y):100j]
    from scipy.interpolate import griddata
    grid_z = griddata((x, y), z, (grid_x, grid_y), method='linear')

    from matplotlib.colors import LogNorm
    print( f"min: {grid_z.min()} max: {grid_z.max()} " )
    contour = ax2.pcolormesh(grid_x, grid_y, grid_z, 
                           norm=LogNorm(vmin=1e-12, vmax=1e-6),
                           cmap='spring', 
                           shading='auto')    

    log_min = 1e-12
    log_max = 1e-6
    levels = 10**np.linspace(-12, -3, 10)
    
    contour_lines = ax2.contour(grid_x, grid_y, grid_z, levels=levels,
                              colors='k', linewidths=0.5)
    ax2.clabel(contour_lines, inline=True, fontsize=8, fmt='%.1e')

    cbar = plt.colorbar(contour)
    ax2.set_xscale('log')

    ax2.set_xlabel(r"$\sigma_{yy}$ (MPa)")
    ax2.set_ylabel("Temperature (°C)")

    for ax in [ax1, ax2] :
        ll = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 20, 25, 30, 40, 55]
        ax.set_xticks(ll)
        ax.set_xticklabels([str(i) for i in ll])
        ax.minorticks_off()

    ax2.set_title("Creep strain rate (1/s)")

    # Plot experimental data
    import itertools
    colors = ['red', 'blue', 'green', 'orange', 'purple', 'yellow', 'pink', 'brown', 'gray', 'cyan', 'magenta', 'lime', 'indigo', 'violet', 'turquoise']
    color_cycle = itertools.cycle(colors)
    df = pd.read_csv("py/csv/fig_5_17_T130.csv", sep="\t")
    l1 = ax1.scatter(df.sig_MPa, df.eps_ss_rate_h/60/60 , c=next(color_cycle), s=15, label=r'Data: $T=130$ °C') 
    df = pd.read_csv("py/csv/fig_5_17_T86.csv", sep="\t")
    l1 = ax1.scatter(df.sig_MPa, df.eps_ss_rate_h/60/60 , c=next(color_cycle), s=15, label=r'Data: $T=86$ °C') 
    df = pd.read_csv("py/csv/fig_5_17_T43.csv", sep="\t")
    l1 = ax1.scatter(df.sig_MPa, df.eps_ss_rate_h/60/60 , c=next(color_cycle), s=15, label=r'Data: $T=43$ °C') 
    ax1.legend()

    fig.tight_layout()

##
#
#
def process_model( rdir, df ) :
    model = f"{rdir}/model/MODEL"
    if not os.path.exists( model ) : 
        wlog(f"    Cannot find model file at '{model}' ... ")

    temperature = None
    sig = None
    sig_x = None
    
    t_pat = re.compile(r'DOMAIN\s+TEST_FRAME\s+T\s+(\d+(?:\.\d+)?)')
    s_pat = re.compile(r'BOUNDARY\s+YP\s+SYY\s+([-+]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)')
    sx_pat = re.compile(r'BOUNDARY\s+XP\s+SXX\s+([-+]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)')
    with open(model, 'r') as file:
        for line in file:
            if m := t_pat.search(line): temperature = float(m.group(1))
            if m := s_pat.search(line): sig = float(m.group(1))
            if m := sx_pat.search(line): sig_x = float(m.group(1))

    df["T"] = temperature
    df["sig"] = sig - sig_x
    print(sig-sig_x)


##
#
#
def process_scalars(rdir) :
    scalars = f"{rdir}/run/csv/scalars.csv"
    if not os.path.exists( scalars ) : 
        wlog(f"    Cannot find scalars file at '{scalars}' ... ")
        return pd.DataFrame()

    df = pd.read_csv(scalars, sep="\t")

    # Select only tri_dy and calculate strain and strain rate
    df = df.pivot_table( values="Value", index=["Time(s)"], columns=['Var'] )
    df=df.reset_index()

    df["eps_yy"] = df.tri_dy / 10  # 10 is the cube size
    df = df[df["Time(s)"]>0]
    x = df["Time(s)"]
    y = df["eps_yy"]
    if not len(x) : return df
    slope, intercept = np.polyfit(x, y, 1)
    df["eps_yy_fit"] = intercept + slope * x
    df["eps_yy_rate"] = slope

    df["rdir"] = rdir

    return df

##
#
#
def process_slurm_out(rdir):
    reg = { 'runtime' : -1 , 'rdir' : rdir }

    slurms = glob(f"{rdir}/slurm*.out")
    if not len(slurms) :
        wlog("    No slurm output found.")
        return reg

    slurm = slurms[-1]
    if len(slurms) != 1 : wlog(f"    Multiple slurm outputs found. Using {slurm} ...")

    n_lines=100
    pattern = re.compile(r'Totals:\s+\d+\s+(\d+\.?\d*)')
    content = subprocess.run(["tail", f"-{n_lines}", slurm], capture_output=True, text=True, check=True).stdout
    m = pattern.findall(content)
    if m: reg["runtime"] = float(m[0])

    return reg


#
#
#
def analytical_model( T_C, sig_MPa ) :

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
        eps_ss_rate += np.exp( -(Q/R/T)**fud  ) * ( sig/sig0ss[i] )**N[i]

    return eps_ss_rate
        
run_main()
plt.show()
