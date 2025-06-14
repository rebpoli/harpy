#!/usr/bin/env -S python 

import os, time

from numpy import log10, logspace, linspace
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import shutil

from env import ilog, flog, elog, dlog, mkdir

all_sig  = - logspace(log10(1e6), log10(50e6), num=15)
all_temp = np.array( [ 46, 86, 126 ] ) + 273  #linspace( 20, 120, num=15 ) + 273

print(all_sig)
print(all_temp)

#
#
#
NEXT_RUN = 0
def create_batch( sig, temp ) :
    global NEXT_RUN
    NEXT_RUN += 1
    run_dir = f"batch/run.{NEXT_RUN}"
    dlog(1 , f"Processing run at '{run_dir}' ...") 

    if os.path.isdir(run_dir): shutil.rmtree( run_dir )

    model_dir = f"{run_dir}/model"
    shutil.copytree( "model_tpl/", model_dir, dirs_exist_ok=True )

    shutil.copy2( "py/run_atena.sh" , run_dir )
    shutil.copy2( "Makefile" , run_dir )

    model_fn = model_dir + "/MODEL"
    
    # Process the template
    subs = { 'TEMPERATURE' : temp,
            'SXX'         : "-10e6",
            'SYY'         : "%.4e"%(sig-10e6),
            'SZZ'         : "-10e6" }
    with open(model_fn) as f: content = f.read()
    for k, v in subs.items():
        content = content.replace(f'%{k}%', str(v))
    with open(model_fn, 'w') as f: f.write(content)
    #

mkdir("batch")
shutil.copy2( "Makefile.batch" , "batch/Makefile" )

for sig in all_sig :
    for temp in all_temp :
        create_batch( sig, temp )
