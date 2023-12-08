#!/bin/bash

#SBATCH --job-name=PIC_OMP1
#SBATCH --time=07:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128

#SBATCH --account=e723-neptune             
#SBATCH --partition=standard
#SBATCH --qos=standard

echo "Start date and time: $(date +"%Y-%m-%d %H:%M:%S")"

export OMP_PLACES=cores
export OMP_PROC_BIND=spread

export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

module load PrgEnv-gnu
module load cray-hdf5-parallel
export PETSC_INSTALL_PATH=/work/e723/e723/csrcnj/lib_install/petsc-3.20.0
export OPPIC_PATH=/work/e723/e723/csrcnj/phd/OP-PIC/lib_oppic
export LD_LIBRARY_PATH=$PETSC_INSTALL_PATH/lib:$LD_LIBRARY_PATH

runFolder=$PWD"/OMP_"$(date +"D_%Y_%m_%d_T_%I_%M_%S")
echo "Creating running folder -> " $runFolder

binPath='/work/e723/e723/csrcnj/phd/OP-PIC/fempic_mpi/bin/'
binName="omp"
echo "Using Binary -> " $binName

hdfOriginalFolder=/work/e723/e723/csrcnj/phd/box_mesh_gen/dats

configFile="/box_archer_omp.param"
file=$PWD'/..'$configFile

# 12000 24000 48000 96000 192000

for run in 1 2; do
    for config in 12000 24000 48000; do
        for tasks in 128 64 32 16 8 4 2 1; do 
            
            export OMP_NUM_THREADS=${tasks}

            folder=$runFolder/$config"_omp"

            echo "Running OMP" $file $config $folder

            mkdir -p $folder
            cp $file $folder
            currentfilename=$folder/$configFile

            copyCommand=$hdfOriginalFolder'/'$config' '$folder'/'
            echo $copyCommand
            cp -r $copyCommand
            escaped_folder="${folder//\//\\/}"
            sed -i "s/STRING global_mesh = \/work\/e723\/e723\/csrcnj\/phd\/box_mesh_gen\/dats\/12000\/mesh.dat/STRING global_mesh = ${escaped_folder}\/${config}\/mesh.dat/" ${currentfilename}
            sed -i "s/STRING inlet_mesh  = \/work\/e723\/e723\/csrcnj\/phd\/box_mesh_gen\/dats\/12000\/inlet.dat/STRING inlet_mesh  = ${escaped_folder}\/${config}\/inlet.dat/" ${currentfilename}
            sed -i "s/STRING wall_mesh   = \/work\/e723\/e723\/csrcnj\/phd\/box_mesh_gen\/dats\/12000\/wall.dat/STRING wall_mesh   = ${escaped_folder}\/${config}\/wall.dat/" ${currentfilename}

            # ---------------------
            command=${binPath}${binName}' '${currentfilename}' >'$folder'/log_N'${tasks}'_C'${config}'_D10_R'${run}'.log 2>&1'
            echo "RUNNING ->"$command
            srun $command | tee $folder/log_N${tasks}_C${config}_D10_R${run}.log;

            sed -i "s/REAL plasma_den     = 1e18/REAL plasma_den     = 0.5e18/" ${currentfilename}
            command=${binPath}${binName}' '${currentfilename}' | tee '$folder'/log_N'${tasks}'_C'${config}'_D5_R'${run}'.log;'
            echo "RUNNING ->"$command
            srun $command | tee $folder/log_N${tasks}_C${config}_D5_R${run}.log;

            sed -i "s/REAL plasma_den     = 0.5e18/REAL plasma_den     = 1e17/" ${currentfilename}
            command=${binPath}${binName}' '${currentfilename}' | tee '$folder'/log_N'${tasks}'_C'${config}'_D1_R'${run}'.log;'
            echo "RUNNING ->"$command
            srun $command | tee $folder/log_N${tasks}_C${config}_D1_R${run}.log;
            # ---------------------

            rmCommand=$folder'/'$config
            echo $rmCommand
            rm -r $rmCommand
        done
    done
done

echo "simulation done"

echo "End date and time: $(date +"%Y-%m-%d %H:%M:%S")"