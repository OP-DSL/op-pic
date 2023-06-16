#!/bin/bash --login
#SBATCH --job-name=SLM_PIC
#SBATCH --partition=v100
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --time=10:00:00

binFolder="/ext-home/zl/phd/OP-PIC/fempic_mpi/bin"
runFolder=$PWD"/LogOMP_"$(date +"D_%Y_%m_%d_T_%I_%M_%S")
echo "creating running folder" $runFolder

source /ext-home/zl/phd/OP-PIC/scripts/source_oneapi

export OMP_PLACES=cores
export OMP_PROC_BIND=close

configs=("coarse_5.param")

for i in ${!configs[@]}; do
    config=${configs[$i]}
    echo $config
    
    file=$PWD/configs/$config
    echo $file

    # ****************************************
    echo "Running OpenMP"

    folder=$runFolder/$config"_omp"

    mkdir -p $folder
    cp $file $folder

    for thr in 48 36 24 12 6 4 2 1; do
        export OMP_NUM_THREADS=${thr}

        echo "OMP " $thr " TEST START"

        srun $binFolder/omp $file | tee $folder/log_thr${thr}.log;	
        # $binFolder/omp $file > $folder/log_thr${thr}.log;

        echo "OMP " $thr " TEST END"
    done
    # ****************************************
done

echo "simulation done"
