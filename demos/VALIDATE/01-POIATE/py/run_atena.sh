#!/bin/bash
#SBATCH --partition cpu   # Project/Partition
#SBATCH -A gger           # Project/Account
#SBATCH -J harpy
#SBATCH --nodes=1         # Number of nodes
#SBATCH --time=1:00:00  # Runtime of this jobs
#SBATCH --cpus-per-task=24


# # #!/bin/bash
# #SBATCH --job-name=poromec
# #SBATCH --nodes=1
# #SBATCH --ntasks=1
# #SBATCH --ntasks-per-node=8
# #SBATCH -A geomec
# #SBATCH --time=12:00:0
# #SBATCH --chdir=/dfs_geral_ep/res/santos/unbs/gger/er/er01/USR/bfq9/lm23_hdd/dev/chimas4d/demos/VALIDATION/19-Monte-Carlo-Constant-Skempton-01/batch

echo START AT
date

# # Load the module environment suitable for the job
# module purge

# pwd

# export SIF=/dfs_geral_ep/res/santos/unbs/gger/er/er01/USR/bfq9/sif/harpy.sif 
# export WORK=/dfs_geral_ep/res/santos/unbs/gger/er/er01/USR/bfq9/work/harpy

# echo SIF:  $SIF
# echo WORK: $WORK
# echo WD:   $WD

# mpirun -np 2 -mca pml ucx -x UCX_TLS=sm,cuda_copy -report-bindings --tag-output \
# singularity exec --nv -B $WORK:/work -B /run/user/$UID --pwd $WD $SIF make run-opt
