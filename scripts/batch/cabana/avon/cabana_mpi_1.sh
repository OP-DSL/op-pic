#!/bin/bash --login
#SBATCH --job-name=PIC_1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3700
#SBATCH --partition=compute

echo "Start date and time: $(date +"%Y-%m-%d %H:%M:%S")"

runFolder=$PWD"/MPI1_"$(date +"D_%Y_%m_%d_T_%I_%M_%S")
echo "Creating running folder" $runFolder

binary="/home/dcs/csrcnj/phd/OP-PIC/cabana_mpi/bin/mpi"
echo "Using Binary" $binary

# module load GCC/10.3.0  OpenMPI/4.1.1

module purge
module load intel/2022a

# export OMP_PLACES=cores
# export OMP_PROC_BIND=close
# export OMP_NUM_THREADS=1

num_nodes=$SLURM_JOB_NUM_NODES

configFile="cabana.param"
file=$PWD/$configFile

for config in 750 1500 3000; do
    for run in 1 2; do

        echo $file $config

        # ****************************************
        echo "Running MPI"

        folder=$runFolder/$config"_mpi"

        mkdir -p $folder
        cp $file $folder
        currentfilename=$folder/$configFile

        sed -i "s/INT nz = 30/INT nz = 60/" ${currentfilename}
        sed -i "s/INT num_part_per_cell = 1500/INT num_part_per_cell = ${config}/" ${currentfilename}
        sed -i "s/STRING cluster = pencil/STRING cluster = block/" ${currentfilename}
        srun $binary $currentfilename  | tee $folder/log_N${num_nodes}_Block_D${config}_R${run}.log;
    done
done

for config in 750 1500 3000; do
    for run in 1 2; do

        echo $file $config

        # ****************************************
        echo "Running MPI"

        folder=$runFolder/$config"_mpi"

        mkdir -p $folder
        cp $file $folder
        currentfilename=$folder/$configFile

        sed -i "s/INT nz = 30/INT nz = 60/" ${currentfilename}
        sed -i "s/INT num_part_per_cell = 1500/INT num_part_per_cell = ${config}/" ${currentfilename}
        srun $binary $currentfilename  | tee $folder/log_N${num_nodes}_Pencil_D${config}_R${run}.log;
    done
done

for config in 750 1500 3000; do
    for run in 1 2; do

        echo $file $config

        # ****************************************
        echo "Running MPI"

        folder=$runFolder/$config"_mpi"

        mkdir -p $folder
        cp $file $folder
        currentfilename=$folder/$configFile

        sed -i "s/INT nz = 30/INT nz = 60/" ${currentfilename}
        sed -i "s/INT num_part_per_cell = 1500/INT num_part_per_cell = ${config}/" ${currentfilename}
        sed -i "s/STRING cluster = pencil/STRING cluster = cart/" ${currentfilename}
        srun $binary $currentfilename  | tee $folder/log_N${num_nodes}_Cart_D${config}_R${run}.log;
    done
done

echo "simulation done"

echo "End date and time: $(date +"%Y-%m-%d %H:%M:%S")"
