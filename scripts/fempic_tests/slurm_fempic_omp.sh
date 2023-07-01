#!/bin/bash --login
#SBATCH --job-name=PIC_OMP
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --time=10:00:00

# binFolder="/ext-home/zl/phd/OP-PIC/fempic_mpi/bin"
binFolder="/home/dcs/csrcnj/phd/OP-PIC/fempic_mpi/bin"
runFolder=$PWD"/LogOMP_"$(date +"D_%Y_%m_%d_T_%I_%M_%S")
echo "creating running folder" $runFolder

# source /ext-home/zl/phd/OP-PIC/scripts/source_oneapi
module load GCC/10.3.0  OpenMPI/4.1.1
module load PETSc/3.15.1

export PETSC_INSTALL_PATH=/scrtp/avon/eb/software/PETSc/3.15.1-foss-2021a
export OPPIC_PATH=/home/dcs/csrcnj/phd/OP-PIC/lib_oppic

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
