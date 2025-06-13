#!/bin/bash
#SBATCH --partition cpu   # Project/Partition
#SBATCH -A gger           # Project/Account
#SBATCH -J harpy
#SBATCH --nodes=1         # Number of nodes
#SBATCH --ntasks=24
#SBATCH --time=1:00:00  # Runtime of this jobs
#SBATCH --cpus-per-task=1

module purge
module load openmpi

echo START AT
date

# # Load the module environment suitable for the job
# module purge

pwd

export SIF=/dfs_geral_ep/res/santos/unbs/gger/er/er01/USR/bfq9/sif/harpy.sif 
export WORK=/dfs_geral_ep/res/santos/unbs/gger/er/er01/USR/bfq9/work/harpy

PWD=$(pwd)
export WD="${PWD/$WORK/\/work}"

echo SIF:  $SIF
echo WORK: $WORK
echo WD:   $WD

## Cleanup
rm -rf run

# mpirun -np $SLURM_NTASKS \
singularity exec --nv -B $WORK:/work --pwd $WD $SIF make _run_
