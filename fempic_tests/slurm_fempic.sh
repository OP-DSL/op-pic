#!/bin/bash --login
#SBATCH --job-name=SLURM_PIC
#SBATCH --partition=v100
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --time=08:00:00

binFolder="/ext-home/zl/phd/OP-PIC/fempic_new/bin1"
runFolder=$PWD"/Log_"$(date +"D_%Y_%m_%d_T_%I_%M_%S")
echo "creating running folder" $runFolder

source /ext-home/zl/phd/OP-PIC/scripts/source_file

configs=("coarse_1.param" "coarse_2.param" "coarse_3.param" "coarse_4.param" "coarse_5.param" 
         "default_1.param" "default_2.param" 
         "coarse_1g.param" "coarse_2g.param" "coarse_3g.param" "coarse_4g.param" "coarse_5g.param" 
         "medium_1.param" "medium_2.param" "medium_3.param" "medium_4.param" "medium_5.param" 
         "medium_1g.param" "medium_2g.param" "medium_3g.param" "medium_4g.param" "medium_5g.param")

for i in ${!configs[@]}; do
    config=${configs[$i]}
    echo $config
    
    file=$PWD/configs/$config
    echo $file

    # ****************************************
    echo "Running Sequential"

    folder=$runFolder/$config"_seq"

    mkdir -p $folder
    cp $file $folder

    export OMP_NUM_THREADS=1
    export OMP_PROC_BIND=close

    srun $binFolder/seq $file | tee $folder/log_seq.log;

    # ****************************************
    echo "Running OpenMP"

    folder=$runFolder/$config"_omp"

    mkdir -p $folder
    cp $file $folder

    for thr in {1,2,4,8,12,16,24,32,40,48}; do
        export OMP_NUM_THREADS=${thr}
        export OMP_PROC_BIND=close

        srun $binFolder/omp $file | tee $folder/log_thr${thr}.log;	

        echo "OMP " $thr " TEST END"
    done
    # ****************************************

done

echo "simulation done"
