#!/usr/bin/env -S python -i

import os, glob
from glob import glob
from os.path import abspath, isdir
import time 
import my_ssh as ssh

start_time = time.time()


rdirs = [ d for d in glob('run.*/') if isdir(d) ]

n_jobs = len(rdirs)
launch_i=0

# Group the targets
N = 1

# Number of launchers
launch_n = int( n_jobs / N )
print(f"Total jobs: {n_jobs}")
print(f"Number of launches: {launch_n}")

# Maximum allocated nodes
MAX = 100

cwd = os.path.abspath(os.getcwd())

# Sessionto launch jobs
ssh.connect( "atenadev01" )


#
# PROCEDURE : Map from container to the outer world (sbatch only works on the outer)
#
import re
from_to = { "^/work/" : "/dfs_geral_ep/res/santos/unbs/gger/er/er01/USR/bfq9/work/harpy/" }
for f,t in from_to.items() :
    cwd = re.sub(f,t,cwd)

all_trg = [];
all_trg.extend( [ (x, "launch-athena") for x in rdirs ]  )
print(f"CWD:{cwd}")

print(f"ALL TARGETS ::: ({len(all_trg)})")
print (all_trg)

# # Clean targets that are up to date (they have a DONE in their slurm file)
# targets = []
# for dname, trg in all_trg :
#     trg_dir = f"{dname}/{trg}"
#     sfiles = glob.glob(f"{trg_dir}/slurm*.out")
#     # Must have a single slurm file, and a DONE in the end
#     if len(sfiles) == 1 : 
#         sfile = sfiles[0]
#         import subprocess
#         cmd = f'tail -n 5 "{sfile}" | grep -q "DONE AT"'
#         result = subprocess.run(cmd, shell=True)
#         if result.returncode == 0:
#             print(f"Skiping {sfile} ... ({result})")
#             continue # Skip file
#     targets.append( ( dname, trg ) )

# print(f"SELECTED TARGETS ::: ({len(targets)})")
# print (targets)

targets = all_trg

#
# PROCEDURE : Spawn the runs
#
LAUNCH = MAX - ssh.n_jobs_running()
n_trgs = len(targets)
i_trgs = 0
while True :
    if i_trgs == len( targets ) : break
    rundir, trg = targets[i_trgs]
    if ( LAUNCH > 0 ) :
        runtime_s = time.time() - start_time
        print(f"[{LAUNCH}] Launching job {i_trgs}/{n_trgs} ({i_trgs/n_trgs*100:.2f} %) - Runtime:{runtime_s/60:.1f} min...")
        cmd = f"cd {cwd}/{rundir} && make {trg}"
        o,e = ssh.cmd(cmd, 0 )
        LAUNCH -= 1
        i_trgs+=1

    else:
        print("# Waiting for cluster space ...")
        time.sleep(15)
        LAUNCH = MAX - ssh.n_jobs_running()

# Sleep a little bit to get all the commands through before the ssh connection dies
print("Last sleep to get all commands through...")
for i in range(10) :
    ssh.cmd("", 0)
    time.sleep(1)

runtime_s = time.time() - start_time
print(f"# Done. - Total runtime: {runtime_s/60:.1f} min")
