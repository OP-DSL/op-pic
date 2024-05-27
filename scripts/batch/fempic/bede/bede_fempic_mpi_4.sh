#!/bin/bash

#SBATCH --job-name=fem4N
#SBATCH --time=02:00:00
#SBATCH --nodes=4

#SBATCH --account=bdyrk17            
#SBATCH --partition=gpu
#SBATCH --gres=gpu:4

echo "Start date and time: $(date +"%Y-%m-%d %H:%M:%S")"

export OMP_NUM_THREADS=1
# export OMP_PLACES=cores
# export OMP_PROC_BIND=close

module load gcc/12.2 openmpi/4.0.5 cuda/12.0.1 openblas/0.3.10

export PETSC_INSTALL_PATH=/users/csrcnl/lib_install/petsc-3.20.1_rel
export OPPIC_PATH=/users/csrcnl/phd/OP-PIC/lib_oppic
export NVCCFLAGS_ADD='-gencode arch=compute_70,code=sm_70'
export LD_LIBRARY_PATH=/users/csrcnl/lib_install/petsc-3.20.1_rel/lib:$LD_LIBRARY_PATH

runFolder=$PWD"/CUDA_4NODE_"$(date +"D_%Y_%m_%d_T_%I_%M_%S")
echo "Creating running folder -> " $runFolder

binary='/users/csrcnl/phd/OP-PIC/app_fempic/bin/cuda_mpi_hdf5'
echo "Using Binary -> " $binary

hdfOriginalFolder=/users/csrcnl/phd/box_mesh_gen/hdf5N
num_nodes=$SLURM_JOB_NUM_NODES

configFile="box_bede.param"
file=$PWD'/'$configFile

for run in 1 2; do
    for config in 384000 768000 1536000; do
            
        folder=$runFolder/$config"_mpi"
	(( totalGPUs=4*$SLURM_JOB_NUM_NODES ))

        echo "Running CUDA-MPI" $file $config $folder

        mkdir -p $folder
        cp $file $folder
        currentfilename=$folder/$configFile

        copyCommand=$hdfOriginalFolder'/box_'$config'.hdf5 '$folder'/'
        echo $copyCommand
        cp $copyCommand
        escaped_folder="${folder//\//\\/}"
        sed -i "s/STRING hdf_filename = \/users\/csrcnl\/phd\/box_mesh_gen\/hdf5\/box_96000.hdf5/STRING hdf_filename = ${escaped_folder}\/box_${config}.hdf5/" ${currentfilename}

        # ---------------------
        echo "RUNNING d 10e17 ->"$totalGPUs $command
        bede-mpirun --bede-par 1ppg ${binary} ${currentfilename} | tee $folder/log_N${num_nodes}_G${totalGPUs}_C${config}_D10_R${run}.log;

        sed -i "s/REAL plasma_den     = 1e18/REAL plasma_den     = 1.3e18/" ${currentfilename}
        echo "RUNNING d 13e17 ->"$totalGPUs $command
        bede-mpirun --bede-par 1ppg ${binary} ${currentfilename} | tee $folder/log_N${num_nodes}_G${totalGPUs}_C${config}_D15_R${run}.log;

        # ---------------------

        rmCommand=$folder'/box_'$config'.hdf5'
        echo $rmCommand
        rm $rmCommand
    done
done

for run in 1 2; do
    for config in 384000 768000 1536000; do
            
        folder=$runFolder/$config"_mpi"
	(( totalGPUs=4*$SLURM_JOB_NUM_NODES ))

        echo "Running MPI" $file $config $folder

        mkdir -p $folder
        cp $file $folder
        currentfilename=$folder/$configFile

        copyCommand=$hdfOriginalFolder'/box_'$config'.hdf5 '$folder'/'
        echo $copyCommand
        cp $copyCommand
        escaped_folder="${folder//\//\\/}"
        sed -i "s/STRING hdf_filename = \/users\/csrcnl\/phd\/box_mesh_gen\/hdf5\/box_96000.hdf5/STRING hdf_filename = ${escaped_folder}\/box_${config}.hdf5/" ${currentfilename}
        sed -i "s/STRING cluster = block/STRING cluster = mpi-block/" ${currentfilename}

        # ---------------------
        echo "RUNNING d 10e17 ->"$totalGPUs $command
        bede-mpirun --bede-par 1ppg ${binary} ${currentfilename} | tee $folder/log_N${num_nodes}_G${totalGPUs}_X${config}_D10_R${run}.log;

        sed -i "s/REAL plasma_den     = 1e18/REAL plasma_den     = 1.3e18/" ${currentfilename}
        echo "RUNNING d 15e17 ->"$totalGPUs $command
        bede-mpirun --bede-par 1ppg ${binary} ${currentfilename} | tee $folder/log_N${num_nodes}_G${totalGPUs}_X${config}_D15_R${run}.log;

        # ---------------------

        rmCommand=$folder'/box_'$config'.hdf5'
        echo $rmCommand
        rm $rmCommand
    done
done

echo "simulation done"

echo "End date and time: $(date +"%Y-%m-%d %H:%M:%S")"
