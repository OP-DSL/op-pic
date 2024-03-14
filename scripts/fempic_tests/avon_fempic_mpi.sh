#!/bin/bash --login
#SBATCH --job-name=PIC_MPI5
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=48
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3700
#SBATCH --partition=compute

binFolder="/home/dcs/csrcnj/phd/OP-PIC/fempic_mpi/bin"
runFolder=$PWD"/LogMPID_"$(date +"D_%Y_%m_%d_T_%I_%M_%S")
echo "creating running folder" $runFolder
echo "bin folder" $binFolder

module load GCC/10.3.0  OpenMPI/4.1.1
module load PETSc/3.15.1
module load CUDA

# export I_MPI_PMI_LIBRARY=/usr/lib/x86_64-linux-gnu/libpmi.so
# export I_MPI_PMI_LIBRARY=/usr/lib64/pmix/lib/libpmi.so

# export OMP_PLACES=cores
# export OMP_PROC_BIND=close
# export OMP_NUM_THREADS=1

configs=("coarse_2.param")

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

    srun -N 2 -n 96 $binFolder/mpi $file  | tee $folder/log_N_2_n_96.log;

    srun -N 2 -n 48 $binFolder/mpi $file  | tee $folder/log_N_2_n_48.log;

    srun -N 1 -n 48 $binFolder/mpi $file  | tee $folder/log_N_1_n_48.log;

    srun -N 1 -n 24 $binFolder/mpi $file  | tee $folder/log_N_1_n_24.log;

    srun -N 1 -n 1 $binFolder/mpi $file  | tee $folder/log_N_1_n_1.log;

    # srun $binFolder/mpi $file > $folder/log_thr${SLURM_PROCID}.log;

    # for thr in 96 72 48 36 24 12 6 4 2 1; do mpirun -np $thr
    #     echo "MPI " $thr " TEST START"

    #     # srun -N 2 -n $thr $binFolder/mpi $file | tee $folder/log_thr${thr}.log;
    #     mpirun -np $thr $binFolder/mpi $file > $folder/log_thr${thr}.log;
        
    #     echo "MPI ALL TESTS END"
    # done
    # ****************************************
done

echo "simulation done"
