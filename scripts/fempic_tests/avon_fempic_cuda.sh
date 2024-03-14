#!/bin/bash --login
#SBATCH --job-name=SLM_PIC
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=3
#SBATCH --cpus-per-task=16
#SBATCH --time=2:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:quadro_rtx_6000:3

binFolder="/home/dcs/csrcnj/phd/OP-PIC/fempic_mpi/bin"
runFolder=$PWD"/LogCUDA_"$(date +"D_%Y_%m_%d_T_%I_%M_%S")
echo "creating running folder" $runFolder

module load GCC/10.3.0  OpenMPI/4.1.1
module load PETSc/3.15.1
module load CUDA

# export OMP_PLACES=cores
# export OMP_PROC_BIND=close

configs=("coarse_2.param")

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

    for t_p_b in 32 64 128 256 512 1024; do

        cp $file $folder
        currentfilename=$folder/$config
        sed -i "s/INT opp_threads_per_block   = 1024/INT opp_threads_per_block   = ${t_p_b}/" ${currentfilename}

        srun -N 2 -n 6 $binFolder/cuda_mpi $currentfilename | tee $folder/log_N_2_t_p_b_${t_p_b}.log;

        srun -N 1 -n 3 $binFolder/cuda_mpi $currentfilename | tee $folder/log_N_1_t_p_b_${t_p_b}.log;	

        # $binFolder/cuda_mpi $currentfilename > $folder/log_alloc${t_p_b}.log;

        echo CUDA_MPI $t_p_b TEST END
    done
done

echo "simulation done"
