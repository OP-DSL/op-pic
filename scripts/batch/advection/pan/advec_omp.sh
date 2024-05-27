#!/bin/bash --login
#SBATCH --job-name=PIC_1
#SBATCH --ntasks-per-node=1
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=64
#SBATCH --partition=mi100

echo "Start date and time: $(date +"%Y-%m-%d %H:%M:%S")"

export OMP_PLACES=cores
export OMP_PROC_BIND=close

runFolder=$PWD"/OMP_ADV_NODE1_S_"$(date +"D_%Y_%m_%d_T_%I_%M_%S")
echo "Creating running folder -> " $runFolder

binPath='/ext-home/zl/phd/neso_test/OP-PIC/advection_mpi/bin/'
binName="omp"
echo "Using Binary -> " $binName

# module purge
# source /opt/intel/oneapi/2021.3/setvars.sh
# export MPICH_CXX=/opt/intel/oneapi/2021.3/compiler/latest/linux/bin/intel64/icpc
# export MPICH_CC=/opt/intel/oneapi/2021.3/compiler/latest/linux/bin/intel64/icc
# export MPICH_F90=/opt/intel/oneapi/2021.3/compiler/latest/linux/bin/intel64/ifort

configFile="config.param"
file=$PWD'/'$configFile

for run in 1; do
    for config in 32 1024; do
        for tasks in 64 32 16 8 4 1; do

            export OMP_NUM_THREADS=${tasks}

            folder=$runFolder/$config"_omp"

            echo "Running OMP" $file $config $folder

            mkdir -p $folder
            cp $file $folder
            currentfilename=$folder/$configFile

            sed -i "s/INT nx = 32/INT nx = ${config}/" ${currentfilename}
            sed -i "s/INT ny = 32/INT ny = ${config}/" ${currentfilename}
            if [[ $config -eq 1024 ]]; then
                sed -i "s/INT n_particles = 6000000/INT n_particles = 100000000/" ${currentfilename}
            fi

            # ---------------------
            command=${binPath}${binName}' '${currentfilename}
            echo "RUNNING dt = 0.001 ->"$command
            # srun -n $tasks $command | tee $folder/log_N${tasks}_C${config}_E3_R${run}.log;
            $command | tee $folder/log_N${tasks}_C${config}_E3_R${run}.log;

            # sed -i "s/REAL dt = 0.001/REAL dt = 0.01/" ${currentfilename}
            # command=${binPath}${binName}' '${currentfilename}
            # echo "RUNNING dt = 0.01 ->"$command
            # # srun -n $tasks $command | tee $folder/log_N${tasks}_C${config}_E2_R${run}.log;
            # $command | tee $folder/log_N${tasks}_C${config}_E2_R${run}.log;

            # sed -i "s/REAL dt = 0.01/REAL dt = 0.1/" ${currentfilename}
            # command=${binPath}${binName}' '${currentfilename}
            # echo "RUNNING dt = 0.1 ->"$command
            # # srun -n $tasks $command | tee $folder/log_N${tasks}_C${config}_E1_R${run}.log;
            # $command | tee $folder/log_N${tasks}_C${config}_E1_R${run}.log;
            # ---------------------
        done
    done
done

echo "simulation done"

echo "End date and time: $(date +"%Y-%m-%d %H:%M:%S")"