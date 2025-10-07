#!/bin/bash
#SBATCH -A gger          
#SBATCH -J harpy_dev
#SBATCH --partition cpu  
#SBATCH --nodes=1       
#SBATCH --time=96:00:00  
#SBATCH --ntasks-per-node=40
#SBATCH --ntasks=40
#SBATCH --exclusive

####SBATCH --gres=gpu:4
####SBATCH --cpus-per-task=10


date
echo PWD:
pwd
echo TC_DIR:
echo $TC_DIR


module purge
module load openmpi


sy() {
    SIF=/dfs_geral_ep/res/santos/unbs/gger/er/er01/USR/bfq9/sif/harpy-2025-09-23.sif

    singularity run \
         -B   /dfs_geral_ep/res/santos/unbs/gger/er/er01/USR/bfq9/work/harpy:/work \
        --pwd $TC_DIR $SIF $@
}



sy make -j20 run
