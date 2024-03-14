#!/bin/bash --login
#SBATCH --job-name=SLM_PIC
#SBATCH --partition=v100
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --time=10:00:00

binFolder="/ext-home/zl/phd/OP-PIC/fempic_mpi/bin"
runFolder=$PWD"/LogCUDA_"$(date +"D_%Y_%m_%d_T_%I_%M_%S")
echo "creating running folder" $runFolder

source /ext-home/zl/phd/OP-PIC/scripts/source_oneapi
module load cuda/toolkit-10.2.89

export OMP_PLACES=cores
export OMP_PROC_BIND=close

configs=("coarse_5.param")

for i in ${!configs[@]}; do
    config=${configs[$i]}
    echo $config
    
    file=$PWD/configs/$config
    echo $file

    # ****************************************
    echo "Running CUDA"

    export OMP_NUM_THREADS=1
    folder=$runFolder/$config"_cuda"
    mkdir -p $folder  

    for t_p_b in 32 64 128 256 512; do
        
        cp $file $folder
        currentfilename=$folder/$config
        sed -i "s/INT opp_threads_per_block   = 256/INT opp_threads_per_block   = ${t_p_b}/" ${currentfilename}

        # srun $binFolder/cuda $currentfilename | tee $folder/log_t_p_b_${t_p_b}.log;	
        $binFolder/cuda $currentfilename > $folder/log_alloc${t_p_b}.log;

        echo CUDA $t_p_b TEST END
    done
done

echo "simulation done"
