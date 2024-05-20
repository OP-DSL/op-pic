#!/bin/bash

#SBATCH --job-name=csrcnj
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1

#SBATCH --account=e723-neptune             
#SBATCH --partition=standard
#SBATCH --qos=standard

echo "Start date and time: $(date +"%Y-%m-%d %H:%M:%S")"

export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=close

export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

module load PrgEnv-gnu
module load cray-hdf5-parallel

runFolder=$PWD"/MPI_"$SLURM_JOB_NUM_NODES"_"$(date +"D_%Y_%m_%d_T_%I_%M_%S")
echo "Creating running folder -> " $runFolder

binary='/mnt/lustre/a2fs-work2/work/e723/e723/csrcnj/phd/OP-PIC/advection_mpi/bin/mpi'
echo "Using Binary -> " $binary

num_nodes=$SLURM_JOB_NUM_NODES

configFile="config.param"
file=$PWD'/'$configFile

for run in 1 2; do
    for config in 128 256 1024 2048; do
           
        folder=$runFolder/$config"_mpi"

        echo "Running MPI" $file $config $folder

        mkdir -p $folder
        cp $file $folder
        currentfilename=$folder/$configFile

        sed -i "s/INT nx = 32/INT nx = ${config}/" ${currentfilename}
        sed -i "s/INT ny = 32/INT ny = ${config}/" ${currentfilename}
        sed -i "s/INT n_particles = 6000000/INT n_particles = 100000000/" ${currentfilename}

        # ---------------------
        srun ${binary} ${currentfilename} | tee $folder/log_N${num_nodes}_C${config}_R${run}.log;
        # ---------------------
    done
done

echo "simulation done"

echo "End date and time: $(date +"%Y-%m-%d %H:%M:%S")"