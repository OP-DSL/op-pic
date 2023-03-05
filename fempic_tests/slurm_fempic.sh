#!/bin/bash --login
#SBATCH --job-name=SLM_PIC
#SBATCH --partition=v100
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --time=10:00:00

binFolder="/ext-home/zl/phd/OP-PIC/fempic_new/bin"
runFolder=$PWD"/Log_"$(date +"D_%Y_%m_%d_T_%I_%M_%S")
echo "creating running folder" $runFolder

source /ext-home/zl/phd/OP-PIC/scripts/source_oneapi
module load cuda/toolkit-10.2.89

export OMP_PLACES=cores
export OMP_PROC_BIND=close

configs=("coarse_1.param" "coarse_2.param" "coarse_3.param" "coarse_4.param" "coarse_5.param" 
         "medium_1.param" "medium_2.param" "medium_3.param" "medium_4.param" "medium_5.param")

for i in ${!configs[@]}; do
    config=${configs[$i]}
    echo $config
    
    file=$PWD/configs/$config
    echo $file

    # ****************************************
    echo "Running CUDA"

    for thr in 1 12 24 48; do

        export OMP_NUM_THREADS=${thr}
        folder=$runFolder/$config"_cuda_"${thr}
        mkdir -p $folder  

        for alloc in 1 5 10 20 40; do
            
            cp $file $folder
            currentfilename=$folder/$config
            sed -i "s/INT opp_allocation_multiple = 1/INT opp_allocation_multiple = ${alloc}/" ${currentfilename}

            srun $binFolder/cuda $currentfilename | tee $folder/log_alloc${alloc}.log;	
            # $binFolder/cuda $currentfilename > $folder/log_alloc${alloc}.log;

            echo CUDA $thr TEST END
        done
    done

    # ****************************************
    echo "Running OpenMP"

    folder=$runFolder/$config"_omp"

    mkdir -p $folder
    cp $file $folder

    for thr in 1 2 4 6 12 24 36 48; do
        export OMP_NUM_THREADS=${thr}

        srun $binFolder/omp $file | tee $folder/log_thr${thr}.log;	
        # $binFolder/omp $file > $folder/log_thr${thr}.log;

        echo "OMP " $thr " TEST END"
    done
    # ****************************************

    echo "Running Sequential"

    folder=$runFolder/$config"_seq"

    mkdir -p $folder
    cp $file $folder

    export OMP_NUM_THREADS=1

    srun $binFolder/seq $file | tee $folder/log_seq.log;
    $binFolder/seq $file > $folder/log_seq.log;

    # ****************************************
done

echo "simulation done"
