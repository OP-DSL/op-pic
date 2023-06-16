#!/bin/bash --login
#SBATCH --job-name=SLM_PIC
#SBATCH --partition=v100
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1

binFolder="/ext-home/zl/phd/OP-PIC/fempic_mpi/bin"
runFolder=$PWD"/LogMPI_"$(date +"D_%Y_%m_%d_T_%I_%M_%S")
echo "creating running folder" $runFolder

source /ext-home/zl/phd/OP-PIC/scripts/source_oneapi

export I_MPI_PMI_LIBRARY=/usr/lib/x86_64-linux-gnu/libpmi.so

export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_NUM_THREADS=1

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

    for thr in 48 36 24 12 6 4 2 1; do mpirun -np $thr
        echo "MPI " $thr " TEST START"

        srun -n $thr $binFolder/mpi $file | tee $folder/log_thr${thr}.log;
        # mpirun -np $thr $binFolder/mpi $file > $folder/log_thr${thr}.log;
        
        echo "MPI ALL TESTS END"
    done
    # ****************************************
done

echo "simulation done"
