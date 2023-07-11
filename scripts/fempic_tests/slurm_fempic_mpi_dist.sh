#!/bin/bash --login
#SBATCH --job-name=PIC_MPID
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=48
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1

thr=768

# binFolder="/ext-home/zl/phd/OP-PIC/fempic_mpi/bin"
binFolder="/home/dcs/csrcnj/phd/OP-PIC/fempic_mpi/bin"
runFolder=$PWD"/LogMPID_"$(date +"D_%Y_%m_%d_T_%I_%M_%S")
echo "creating running folder" $runFolder
echo "bin folder" $binFolder

# source /ext-home/zl/phd/OP-PIC/scripts/source_oneapi

module load GCC/10.3.0  OpenMPI/4.1.1
module load PETSc/3.15.1
module load pmi

export PETSC_INSTALL_PATH=/scrtp/avon/eb/software/PETSc/3.15.1-foss-2021a
export OPPIC_PATH=/home/dcs/csrcnj/phd/OP-PIC/lib_oppic

# export I_MPI_PMI_LIBRARY=/usr/lib/x86_64-linux-gnu/libpmi.so
export I_MPI_PMI_LIBRARY=/usr/lib64/pmix/lib/libpmi.so

# export OMP_PLACES=cores
# export OMP_PROC_BIND=close
# export OMP_NUM_THREADS=1

configs=("coarse_5.param")

for i in ${!configs[@]}; do
    config=${configs[$i]}
    echo $config
    
    file=$PWD/configs/$config
    echo $file

    # ****************************************
    echo "Running MPI"

    folder=$runFolder/$config"_mpi"

    mkdir -p $folder
    cp $file $folder

    srun $binFolder/mpi $file | tee $folder/log_thr${thr}.log;

    # for thr in 96 72 48 36 24 12 6 4 2 1; do mpirun -np $thr
    #     echo "MPI " $thr " TEST START"

    #     # srun -N 2 -n $thr $binFolder/mpi $file | tee $folder/log_thr${thr}.log;
    #     mpirun -np $thr $binFolder/mpi $file > $folder/log_thr${thr}.log;
        
    #     echo "MPI ALL TESTS END"
    # done
    # ****************************************
done

echo "simulation done"
