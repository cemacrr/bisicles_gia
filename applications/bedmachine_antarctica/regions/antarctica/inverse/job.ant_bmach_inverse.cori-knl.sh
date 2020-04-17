#!/bin/bash 
#
# job script for nersc cori
#
#SBATCH -A m1041
#SBATCH -J ant_bdmach_inverse
##SBATCH --qos=debug
#SBATCH --qos=regular
#SBATCH --time=10:00:00
#SBATCH --nodes=4
#SBATCH --tasks-per-node=68
#SBATCH --constraint=knl


NCORE=$((68 * $SLURM_JOB_NUM_NODES))
echo "n_node: $SLURM_JOB_NUM_NODES , n_core: $NCORE"
export DRIVER=$HOME/cori-bisicles/BISICLES/code/exec2D/driver2d.Linux.64.CC.ftn.DEBUG.OPT.MPI.PETSC.ex

cd $SLURM_SUBMIT_DIR

INFILE=$SLURM_SUBMIT_DIR/inputs.ant.inverse

export CH_TIMER=1
export CH_OUTPUT_INTERVAL=999
export PYTHONPATH=$PWD:$PYTHONPATH

echo "srun $DRIVER $INFILE"
srun $DRIVER $INFILE
