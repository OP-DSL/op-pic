#!/bin/bash --login
#SBATCH --job-name=PIC_1k_sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3700
#SBATCH --partition=compute

echo "Start date and time: $(date +"%Y-%m-%d %H:%M:%S")"

# export OMP_PLACES=cores
# export OMP_PROC_BIND=close

runFolder=$PWD"/MPI_ADV_NODE1_L_"$(date +"D_%Y_%m_%d_T_%I_%M_%S")
echo "Creating running folder -> " $runFolder

binPath='/home/dcs/csrcnj/phd/OP-PIC/advection_mpi/bin/'
binName="mpi"
echo "Using Binary -> " $binName

module purge
module load intel/2022a
export MPICH_CXX=/scrtp/avon/eb/software/intel-compilers/2022.1.0/compiler/2022.1.0/linux/bin/intel64/icpc
export MPICH_CC=/scrtp/avon/eb/software/intel-compilers/2022.1.0/compiler/2022.1.0/linux/bin/intel64/icc
export MPICH_F90=/scrtp/avon/eb/software/intel-compilers/2022.1.0/compiler/2022.1.0/linux/bin/intel64/ifort

CPU_BIND="map_cpu:0,4,8,12,16,20,24,28,32,36,40,44,1,5,9,13,17,21,23,25,29,33,41,45,2,6,10,14,18,22,26,30,34,38,42,46,3,7,11,15,19,27,31,35,37,39,43,47"

configFile="config.param"
file=$PWD'/'$configFile

for run in 1 2; do
    for config in 1024; do
        for tasks in 48 24 12 8 4 1; do

            folder=$runFolder/$config"_mpi"

            echo "Running MPI" $file $config $folder

            mkdir -p $folder
            cp $file $folder
            currentfilename=$folder/$configFile

            sed -i "s/INT nx = 32/INT nx = ${config}/" ${currentfilename}
            sed -i "s/INT ny = 32/INT ny = ${config}/" ${currentfilename}
            sed -i "s/INT n_particles = 6000000/INT n_particles = 100000000/" ${currentfilename}

            # ---------------------
            command=${binPath}${binName}' '${currentfilename}
            echo "RUNNING dt = 0.001 ->"$command
            srun -n $tasks --cpu-bind=${CPU_BIND} $command | tee $folder/log_N${tasks}_C${config}_E3_R${run}.log;

            sed -i "s/REAL dt = 0.001/REAL dt = 0.01/" ${currentfilename}
            command=${binPath}${binName}' '${currentfilename}
            echo "RUNNING dt = 0.01 ->"$command
            srun -n $tasks --cpu-bind=${CPU_BIND} $command | tee $folder/log_N${tasks}_C${config}_E2_R${run}.log;

            sed -i "s/REAL dt = 0.01/REAL dt = 0.1/" ${currentfilename}
            command=${binPath}${binName}' '${currentfilename}
            echo "RUNNING dt = 0.1 ->"$command
            srun -n $tasks --cpu-bind=${CPU_BIND} $command | tee $folder/log_N${tasks}_C${config}_E1_R${run}.log;
            # ---------------------
        done
    done
done

echo "simulation done"

echo "End date and time: $(date +"%Y-%m-%d %H:%M:%S")"