#!/bin/bash -l
#SBATCH -p debug
#SBATCH -N 1
#SBATCH -t 00:05:00
#SBATCH -J s4cmbrocks
#SBATCH -C haswell

source $HOME/.bashrc.ext
cd $SLURM_SUBMIT_DIR

## Update the path!
path_to_scripts=$PWD

time srun -n 12 python-mpi ${path_to_scripts}/example/test/simple_app.py \
    -inifile ${path_to_scripts}/examples/inifiles/simple_parameters.py -tag run_0
