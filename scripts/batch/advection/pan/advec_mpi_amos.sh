#!/bin/bash --login

#SBATCH --job-name=PIC_1
#SBATCH --ntasks-per-node=48
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1
#SBATCH --partition=mi100

echo "Start date and time: $(date +"%Y-%m-%d %H:%M:%S")"

export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_NUM_THREADS=1

runFolder=$PWD"/AMOS_MPI_ADV_"$(date +"D_%Y_%m_%d_T_%I_%M_%S")
echo "Creating running folder -> " $runFolder

binary='/ext-home/zl/phd/neso_test/OP-PIC/advection_mpi/bin/mpi'
echo "Using Binary -> " $binary

module purge
source /opt/intel/oneapi/2021.3/setvars.sh
# export MPICH_CXX=/opt/intel/oneapi/2021.3/compiler/latest/linux/bin/intel64/icpc
# export MPICH_CC=/opt/intel/oneapi/2021.3/compiler/latest/linux/bin/intel64/icc
# export MPICH_F90=/opt/intel/oneapi/2021.3/compiler/latest/linux/bin/intel64/ifort

configFile="config.param"
file=$PWD'/'$configFile

ppc=1000

for run in 1 2; do
    for config in 32 128 1024; do
        for ppc in 1000 5000 10000 15000; do
            for tasks in 64; do

                folder=$runFolder/$config"_mpi"

                echo "Running MPI" $file $config $folder

                mkdir -p $folder
                cp $file $folder
                currentfilename=$folder/$configFile

                parts=$((config*config*ppc))
                sed -i "s/INT n_particles = 6000000/INT n_particles = ${parts}/" ${currentfilename}

                config=$((config*2))
                sed -i "s/INT nx = 32/INT nx = ${config}/" ${currentfilename}
                sed -i "s/INT ny = 32/INT ny = ${config}/" ${currentfilename}
                
                # ---------------------
                echo "RUNNING dt = 0.001 ->"
                # srun -n $tasks $binary ${currentfilename} | tee $folder/log_N${tasks}_C${config}_E3_R${run}.log;
                mpirun -np $tasks $binary ${currentfilename} | tee $folder/log_N${tasks}_C${config}_P${ppc}_E3_R${run}.log;

                config=$((config/2))
                sed -i "s/REAL dt = 0.001/REAL dt = 0.01/" ${currentfilename}
                echo "RUNNING dt = 0.01 ->"
                # srun -n $tasks $binary ${currentfilename} | tee $folder/log_N${tasks}_C${config}_E2_R${run}.log;
                mpirun -np $tasks $binary ${currentfilename} | tee $folder/log_N${tasks}_C${config}_E2_R${run}.log;

                sed -i "s/REAL dt = 0.01/REAL dt = 0.1/" ${currentfilename}
                echo "RUNNING dt = 0.1 ->"
                # srun -n $tasks $binary ${currentfilename} | tee $folder/log_N${tasks}_C${config}_E1_R${run}.log;
                mpirun -np $tasks $binary ${currentfilename} | tee $folder/log_N${tasks}_C${config}_E1_R${run}.log;
                # ---------------------
            done
        done
    done
done

echo "simulation done"

echo "End date and time: $(date +"%Y-%m-%d %H:%M:%S")"


