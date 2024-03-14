#!/bin/bash --login
#SBATCH --job-name=PIC_OMP
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=3700
#SBATCH --partition=compute

runFolder=$PWD"/OMP_"$(date +"D_%Y_%m_%d_T_%I_%M_%S")
echo "Creating running folder" $runFolder

binary="/home/dcs/csrcnj/phd/OP-PIC/fempic_mpi/bin/omp"
echo "Using Binary" $binary

export OMP_PLACES=cores
export OMP_PROC_BIND=spread
# export OMP_PROC_BIND=close

num_nodes=$SLURM_JOB_NUM_NODES

configFile="box_avon.param"
file=$PWD/../$configFile

# 24000 12000 48000 96000

for run in 1 2; do
    for config in 192000; do
        for omp in 48 36 24 16 8 4 2 1; do

            export OMP_NUM_THREADS=$omp

            echo $file $config

            # ****************************************
            echo "Running OMP"

            folder=$runFolder/$config"_omp"

            mkdir -p $folder
            cp $file $folder
            currentfilename=$folder/$configFile

            sed -i "s/STRING global_mesh = \/home\/dcs\/csrcnj\/phd\/box_mesh_gen\/12000\/mesh.dat/STRING global_mesh = \/home\/dcs\/csrcnj\/phd\/box_mesh_gen\/${config}\/mesh.dat/" ${currentfilename}
            sed -i "s/STRING inlet_mesh = \/home\/dcs\/csrcnj\/phd\/box_mesh_gen\/12000\/inlet.dat/STRING inlet_mesh = \/home\/dcs\/csrcnj\/phd\/box_mesh_gen\/${config}\/inlet.dat/" ${currentfilename}
            sed -i "s/STRING wall_mesh = \/home\/dcs\/csrcnj\/phd\/box_mesh_gen\/12000\/wall.dat/STRING wall_mesh = \/home\/dcs\/csrcnj\/phd\/box_mesh_gen\/${config}\/wall.dat/" ${currentfilename}

            srun $binary $currentfilename  | tee $folder/log_O${omp}_C${config}_D10_R${run}.log;

            sed -i "s/REAL plasma_den     = 1e18/REAL plasma_den     = 0.5e18/" ${currentfilename}

            srun $binary $currentfilename  | tee $folder/log_O${omp}_C${config}_D5_R${run}.log;

            sed -i "s/REAL plasma_den     = 0.5e18/REAL plasma_den     = 1e17/" ${currentfilename}

            srun $binary $currentfilename  | tee $folder/log_O${omp}_C${config}_D1_R${run}.log;
            
        done
    done
done

echo "simulation done"
