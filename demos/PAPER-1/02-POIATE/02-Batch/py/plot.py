#!/usr/bin/env -S python -i

import sys, os
sys.path.append('../py')
sys.stdout.reconfigure(line_buffering=True)

import pandas as pd
import numpy as np

import matplotlib
import matplotlib.pyplot as plt
# matplotlib.use('Agg')  ## Noninteractive

from plot_util import savefig

cwd = os.path.dirname(__file__)

fn = f"{cwd}/my.mplstyle"
plt.style.use(fn)

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
        

df = pd.read_pickle(f'{cwd}/../cache/batch_db.pkl')

SEL_TEMP = np.array([ 130, 86, 46 ]) + 273

import itertools
colors = ['red', 'blue', 'green', 'orange', 'purple', 'yellow', 'pink', 'brown', 'gray', 'cyan', 'magenta', 'lime', 'indigo', 'violet', 'turquoise']
color_cycle = itertools.cycle(colors)

all_axes = []

#
# PLOT 1 - ANALYTICAL VS NUMERICAL
#

fig, ax = plt.subplots(1, 1, figsize=(7,5))
all_axes.append(ax)
groups = list(df.groupby("T", sort=True))
for k, g in reversed(groups):
    if not k in SEL_TEMP : continue
    print(k) 
    if len(SEL_TEMP) and not k in SEL_TEMP : continue
    ax.scatter( -g.sig/1e6, -g.eps_yy_rate, marker='x', lw=1, label=f"Numerical: T={k-273:.1f} °C", c=next(color_cycle))

an_data = { 'T' : [], 'sig' : [], 'eps_rate' : [] }
color_cycle = itertools.cycle(colors)
all_temp = df["T"].drop_duplicates()
for temp in reversed(list(all_temp)) :
    if len(SEL_TEMP) and not temp in SEL_TEMP : continue
    all_rate = []
    all_sig = np.logspace(np.log10(1e6), np.log10(50e6), num=100)
    for sig in all_sig :
        an = analytical_model( temp-273, sig/1e6 )
        all_rate.append( an )
    ax.plot( all_sig/1e6, all_rate, '-', c=next(color_cycle) , label=f"Analytical: $T={temp-273}$ °C")#, c=next(color_cycle) )

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylabel(r"Strain rate (1/s)")
ax.set_xlabel(r"Deviatoric Stress (MPa)")
ax.legend()

savefig(fig, "png/analitycal_vs_numerical.png")

#
# PLOT 2 - NUMERICAL vs EXPERIMENTAL DATA
#

fig, ax = plt.subplots(1, 1, figsize=(7,5))
color_cycle = itertools.cycle(colors)
all_axes.append(ax)
groups = list(df.groupby("T", sort=True))
for k, g in reversed(groups):
    if not k in SEL_TEMP : continue
    print(k) 
    if len(SEL_TEMP) and not k in SEL_TEMP : continue
    ax.plot( -g.sig/1e6, -g.eps_yy_rate, ls='--', lw=1, label=f"Numerical: T={k-273:.1f} °C", c=next(color_cycle))

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylabel(r"Strain rate (1/s)")
ax.set_xlabel(r"Deviatoric Stress (MPa)")


# Plot experimental data
import itertools
colors = ['red', 'blue', 'green', 'orange', 'purple', 'yellow', 'pink', 'brown', 'gray', 'cyan', 'magenta', 'lime', 'indigo', 'violet', 'turquoise']
color_cycle = itertools.cycle(colors)
df_ = pd.read_csv(f"{cwd}/../raw_data/fig_5_17_T130.csv", sep="\t")
l1 = ax.scatter(df_.sig_MPa, df_.eps_ss_rate_h/60/60 , c=next(color_cycle), s=15, label=r'Experimental: $T=130$ °C') 
df_ = pd.read_csv(f"{cwd}/../raw_data/fig_5_17_T86.csv", sep="\t")
l1 = ax.scatter(df_.sig_MPa, df_.eps_ss_rate_h/60/60 , c=next(color_cycle), s=15, label=r'Experimental: $T=86$ °C') 
df_ = pd.read_csv(f"{cwd}/../raw_data/fig_5_17_T43.csv", sep="\t")
l1 = ax.scatter(df_.sig_MPa, df_.eps_ss_rate_h/60/60 , c=next(color_cycle), s=15, label=r'Experimental: $T=43$ °C') 

ax.legend()

savefig(fig, "png/numerical_vs_experimental.png")

#
# PLOT 3 - colormap
#

x = -df["sig"].values / 1e6  # MPa
y = df["T"].values - 273 # ºC
z = -df['eps_yy_rate'].values

an_data = { 'T' : [], 'sig' : [], 'eps_rate' : [] }
color_cycle = itertools.cycle(colors)
full_temp_range = np.arange(10,140,1)+273
for temp in reversed(list(full_temp_range)) :
    all_rate = []
    all_sig = np.logspace(np.log10(1e6), np.log10(50e6), num=100)
    for sig in all_sig :
        an_data['T'].append(temp-273)
        an_data['sig'].append( sig/1e6 )
        an = analytical_model( temp-273, sig/1e6 )
        an_data['eps_rate'].append( an )

# cria a malha
x = an_data['sig']
y = an_data['T']
z = an_data['eps_rate']
            
grid_x, grid_y = np.mgrid[min(x):max(x):100j, min(y):max(y):100j]
from scipy.interpolate import griddata
grid_z = griddata((x, y), z, (grid_x, grid_y), method='linear')

fig, ax2 = plt.subplots(1, 1, figsize=(7,5))
all_axes.append(ax)
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

ax2.set_xlabel(r"Deviatoric stress (MPa)")
ax2.set_ylabel("Temperature (°C)")

for ax in all_axes :
    ll = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 20, 25, 30, 40, 55]
    ax.set_xticks(ll)
    ax.set_xticklabels([str(i) for i in ll])
    ax.minorticks_off()

savefig(fig, "png/colormap.png")

plt.show()
