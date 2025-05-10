#!/usr/bin/env -S python -i

import pandas as pd
import matplotlib.pyplot as plt
import os

fn = "py/paper.mplstyle"
if os.path.isfile(fn) :
    plt.style.use('default')   ## reset!
    plt.style.use(fn)

df = pd.read_csv("run/csv/lin_0.csv", sep="\t")

# Index(['Timestep', 'Time(s)', 'Time(day)', 'Var', 'X', 'Y', 'Z', 'Value'], dtype='object')

fig, ax = plt.subplots()

sxxdf = df[ df.Var == "STOTXX" ]
sxxdf = sxxdf[ sxxdf.X > -0.02 ]
sxxdf = sxxdf[ sxxdf.X <  0 ]
ax.plot( sxxdf["Time(day)"], sxxdf["Value"] )
ax.set_xscale('log')

plt.show()
