#!/bin/bash

#SBATCH --job-name=csrcnj
#SBATCH --time=0:20:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1

#SBATCH --account=e723-neptune             
#SBATCH --partition=standard
#SBATCH --qos=standard

# Set the number of threads to 1
#   This prevents any threaded system libraries from automatically 
#   using threading.
# export OMP_NUM_THREADS=1
# export OMP_PLACES=cores
# export OMP_PROC_BIND=close

# Propagate the cpus-per-task setting from script to srun commands
#    By default, Slurm does not propagate this setting from the sbatch
#    options to srun commands in the job script. If this is not done,
#    process/thread pinning may be incorrect leading to poor performance
export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

# Launch the parallel job
#   Using 512 MPI processes and 128 MPI processes per node
#   srun picks up the distribution from the sbatch options

module load PrgEnv-gnu

module load cray-hdf5-parallel
export PETSC_INSTALL_PATH=/work/e723/e723/csrcnj/lib_install/petsc-3.20.0
export OPPIC_PATH=/work/e723/e723/csrcnj/phd/OP-PIC/lib_oppic

export LD_LIBRARY_PATH=$PETSC_INSTALL_PATH/lib:$LD_LIBRARY_PATH

module load valgrind4hpc/2.12.10
binName='mpi_hdf5'
command='/work/e723/e723/csrcnj/phd/OP-PIC/fempic_mpi/bin/'${binName}' /work/e723/e723/csrcnj/phd/OP-PIC/scripts/fempic_tests/configs/coarse_2.param'

# srun --distribution=block:block --hint=nomultithread /work/e723/e723/csrcnj/phd/OP-PIC/fempic_mpi/bin/mpi_hdf5 /work/e723/e723/csrcnj/phd/OP-PIC/scripts/fempic_tests/configs/coarse_2.param

echo $command

srun $command