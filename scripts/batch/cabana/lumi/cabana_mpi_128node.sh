#!/bin/bash

#SBATCH --job-name=cab_128            # Job name
#SBATCH --output=cab_128.o%j          # Name of stdout output file
#SBATCH --error=cab_128.e%j           # Name of stderr error file
#SBATCH --partition=standard-g    # Partition (queue) name
#SBATCH --nodes=128                   # Total number of nodes 
#SBATCH --ntasks-per-node=8          # 8 MPI ranks per node, 16 total (2x8)
#SBATCH --gpus-per-node=8            # Allocate one gpu per MPI rank
#SBATCH --time=0-00:50:00            # Run time (d-hh:mm:ss)
##SBATCH --mail-type=all             # Send email at begin and end of job
#SBATCH --account=project_465001068  # Project for billing
##SBATCH --mail-user=username@domain.com

echo "Start date and time: $(date +"%Y-%m-%d %H:%M:%S")"

runFolder=$PWD"/MPI_"${SLURM_JOB_NUM_NODES}"_"$(date +"D_%Y_%m_%d_T_%I_%M_%S")
echo "Creating running folder" $runFolder

runFolder=$PWD"/MPI_N"${SLURM_JOB_NUM_NODES}"_SR"${use_seg_red}"_"$(date +"D_%Y_%m_%d_T_%I_%M_%S")
echo "Creating running folder -> " $runFolder

binpath=/users/lantraza/phd/OP-PIC/app_cabanapic/bin
binary=$binpath"/hip_mpi"
echo "Using Binary -> " $binary

configFile="cabana.param"
file=/users/lantraza/phd/OP-PIC/scripts/batch/cabana/lumi/$configFile
echo "Config file -> " $file

cd $binpath
echo "********************************************************"
gitbranch=$(git branch | sed -n -e 's/^\* \(.*\)/\1/p')
gitcommit=$(git log -n 1 $gitbranch)
echo "Git branch " $gitbranch
echo "Git commit " $gitcommit
echo "********************************************************"
cd -

export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=close

CPU_BIND="map_cpu:49,57,17,25,1,9,33,41"
export MPICH_GPU_SUPPORT_ENABLED=1

module load LUMI/23.09
module load partition/G
module load cray-hdf5-parallel
module load SCOTCH/6.1.3-cpeCray-23.09

num_nodes=$SLURM_JOB_NUM_NODES

gpus=8

for config in 750 1500 3000; do
    for run in 1 2; do

        (( totalGPUs=$gpus*$SLURM_JOB_NUM_NODES ))
        echo $file $config

        # ****************************************
        echo "Running MPI BLOCK " $config $run $totalGPUs

        folder=$runFolder/$config"_mpi"

        mkdir -p $folder
        cp $file $folder
        currentfilename=$folder/$configFile
        
        (( rnz=$totalGPUs*60 ))
        sed -i "s/INT nz = 30/INT nz = ${rnz}/" ${currentfilename}
        sed -i "s/INT num_part_per_cell = 1500/INT num_part_per_cell = ${config}/" ${currentfilename}
        sed -i "s/STRING cluster = pencil/STRING cluster = block/" ${currentfilename}
        srun --cpu-bind=${CPU_BIND} -n $totalGPUs $binary $currentfilename  | tee $folder/log_N${num_nodes}_G${totalGPUs}_Block_GD_0_D${config}_R${run}.log;
    done
done

# for config in 750 1500; do
#     for run in 1 2; do

#         echo $file $config

#         # ****************************************
#         echo "Running CUDA-MPI"

#         folder=$runFolder/$config"_mpi"

#         mkdir -p $folder
#         cp $file $folder
#         currentfilename=$folder/$configFile

#         sed -i "s/INT nz = 30/INT nz = 7680/" ${currentfilename}
#         sed -i "s/INT num_part_per_cell = 1500/INT num_part_per_cell = ${config}/" ${currentfilename}
#         sed -i "s/STRING cluster = pencil/STRING cluster = cart/" ${currentfilename}
#         bede-mpirun --bede-par 1ppg $binary $currentfilename  | tee $folder/log_N${num_nodes}_Cart_D${config}_R${run}.log;
#     done
# done

# for config in 750 1500; do
#     for run in 1 2; do

#         echo $file $config

#         # ****************************************
#         echo "Running CUDA-MPI"

#         folder=$runFolder/$config"_mpi"

#         mkdir -p $folder
#         cp $file $folder
#         currentfilename=$folder/$configFile

#         sed -i "s/INT nz = 30/INT nz = 7680/" ${currentfilename}
#         sed -i "s/INT num_part_per_cell = 1500/INT num_part_per_cell = ${config}/" ${currentfilename}
#         bede-mpirun --bede-par 1ppg $binary $currentfilename  | tee $folder/log_N${num_nodes}_Pencil_D${config}_R${run}.log;
#     done
# done


echo "simulation done"

echo "End date and time: $(date +"%Y-%m-%d %H:%M:%S")"
