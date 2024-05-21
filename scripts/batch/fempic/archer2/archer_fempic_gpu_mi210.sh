#!/bin/bash

#SBATCH --job-name=test
#SBATCH --gpus=1
#SBATCH --time=00:10:00
#SBATCH --account=e723-neptune 
#SBATCH --partition=gpu
#SBATCH --qos=gpu-shd

module purge
module load PrgEnv-amd
module load rocm
module load craype-accel-amd-gfx90a
module load craype-x86-milan

export PETSC_INSTALL_PATH=/work/e723/e723/csrcnj/lib_install/petsc-3.20.0_Milan
export LD_LIBRARY_PATH=$PETSC_INSTALL_PATH/lib:$LD_LIBRARY_PATH

srun --ntasks=1 rocm-smi

# srun --ntasks=2 --cpus-per-task=1 /work/e723/e723/csrcnj/phd/OP-PIC/advection_mpi/bin/hip_mpi /work/e723/e723/csrcnj/phd/OP-PIC/advection_mpi/configs/advec.param
srun --ntasks=1 --cpus-per-task=1 /work/e723/e723/csrcnj/phd/OP-PIC/app_fempic/bin/hip_mpi /work/e723/e723/csrcnj/phd/OP-PIC/app_fempic/configs/coarse.param

srun --ntasks=1 rocm-smi
