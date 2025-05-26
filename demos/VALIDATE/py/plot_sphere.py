#!/usr/bin/env -S python

import pandas as pd
import matplotlib.pyplot as plt
import os, time
from matplotlib.animation import FuncAnimation

fn = "py/paper.mplstyle"
if os.path.isfile(fn) :
    plt.style.use('default')   ## reset!
    plt.style.use(fn)


#
#
#

fig, ax = plt.subplots(figsize=(7, 4))

ss_csv = "03-SPHERE-SS/run/csv/lin_0.csv"
ss_df = pd.read_csv(ss_csv, sep="\t")
ss_df = ss_df[ ss_df["Time(s)"]>0 ]
ss_df = ss_df[ ss_df.X > -0.02 ]
ss_df = ss_df[ ss_df.X <  0 ]

sxxdf = ss_df[ ss_df.Var == "STOTXX" ]
syydf = ss_df[ ss_df.Var == "STOTYY" ]
szzdf = ss_df[ ss_df.Var == "STOTZZ" ]
ax.plot( sxxdf["Time(day)"], sxxdf["Value"], ls='-', label=r"$\sigma_{xx}$ (ss)" , c='b')
ax.plot( syydf["Time(day)"], syydf["Value"], ls='-', label=r"$\sigma_{yy}$ (ss)" , c='r')
ax.plot( szzdf["Time(day)"], szzdf["Value"], ls='-', label=r"$\sigma_{zz}$ (ss)" , c='g')

tr_csv = "02-SPHERE/run/csv/lin_0.csv"
tr_df = pd.read_csv(tr_csv, sep="\t")
tr_df = tr_df[ tr_df["Time(s)"]>0 ]
tr_df = tr_df[ tr_df.X > -0.02 ]
tr_df = tr_df[ tr_df.X <  0 ]

sxxdf = tr_df[ tr_df.Var == "STOTXX" ]
syydf = tr_df[ tr_df.Var == "STOTYY" ]
szzdf = tr_df[ tr_df.Var == "STOTZZ" ]
ax.plot( sxxdf["Time(day)"], sxxdf["Value"], ls='--', label=r"$\sigma_{xx}$ (tr)", c='b')
ax.plot( syydf["Time(day)"], syydf["Value"], ls='--', label=r"$\sigma_{yy}$ (tr)", c='r')
ax.plot( szzdf["Time(day)"], szzdf["Value"], ls='--', label=r"$\sigma_{zz}$ (tr)", c='g')


#
#
#
ax.set_ylim( 1e7, 4e7 )
ax.set_xscale('log')
ax.set_xlabel("Time (days)")
ax.set_ylabel(r"Stress at the sphere center (MPa)")
ax.legend()

plt.show()
