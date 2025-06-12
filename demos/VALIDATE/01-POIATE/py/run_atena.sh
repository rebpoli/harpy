#!/bin/bash
#SBATCH --partition gpu   # Project/Partition
#SBATCH -A gger           # Project/Account
#SBATCH -J lbpm_dp
#SBATCH --nodes=1         # Number of nodes
#SBATCH --time=10:00:00  # Runtime of this jobs
# Number of tasks per node
#SBATCH --ntasks-per-node=8
# Number of cpus per task
# Number of gpus
#SBATCH --gres=gpu:4

####SBATCH --cpus-per-task=10

#
# This must be defined in the environment. We may want to add validation or fallback here
#
# export SIF=/dfs_geral_ep/res/santos/unbs/gger/er/er01/USR/bfq9/sif/lbpm.sif
# export TOP_DIR=$PWD/../..

echo START AT
date
export GPU_ARCH="sm70" # Valid options are: sm60, sm70, sm75, sm80

# Load the module environment suitable for the job
module purge
module load openmpi/4.1.4_ucx1.13.1_cuda_11.3_gcc

pwd


echo RUNNING LBPM_MODE=$LBPM_MODE

if [ "$LBPM_MODE" == "singlephase" ] ; then
  echo Running SINGLEPHASE routine.
  cat input_single.db
  mpirun -np 2 -mca pml ucx -x UCX_TLS=sm,cuda_copy -report-bindings --tag-output \
  singularity run --nv -B $TOPDIR $SIF $LBPM_DIR/lbpm_permeability_simulator input_single.db
  cat Permeability.csv
fi

if [ "$LBPM_MODE" == "morphdrain" ] ; then
  # Run morphdrain
  cat input_morphdrain.db
##  mpirun -np 4 -mca pml ucx -x UCX_TLS=sm,cuda_copy -report-bindings --tag-output \
  singularity run --nv -B $TOPDIR $SIF $LBPM_DIR/lbpm_morphdrain_pp input_morphdrain.db
fi

if [ "$LBPM_MODE" == "multiphase" ] ; then
  echo Running MULTIPHASE routine.
  # Run color solver
  cat input_color.db
  mpirun -np 2 -mca pml ucx -x UCX_TLS=sm,cuda_copy -report-bindings --tag-output \
  singularity run --nv -B $TOPDIR $SIF $LBPM_DIR/lbpm_color_simulator input_color.db
fi

if [ "$LBPM_MODE" == "multiphase-restart" ] ; then
  echo Running MULTIPHASE-RESTART routine.
  # Run color solver
  cat Restart.db
  mpirun -np 2 -mca pml ucx -x UCX_TLS=sm,cuda_copy -report-bindings --tag-output \
  singularity run --nv -B $TOPDIR $SIF $LBPM_DIR/lbpm_color_simulator Restart.db 
fi

echo DONE AT
date
#cat ./relperm.csv

