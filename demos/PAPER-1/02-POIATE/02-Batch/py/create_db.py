#!/usr/bin/env -S python -i
import sys, os
sys.path.append('../py')
sys.stdout.reconfigure(line_buffering=True)

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

    df_.to_pickle("cache/batch_db.pkl")


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


run_main()
# plt.show()
