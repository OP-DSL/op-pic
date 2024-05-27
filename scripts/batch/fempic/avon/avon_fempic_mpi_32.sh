#!/bin/bash --login
#SBATCH --job-name=PIC_32
#SBATCH --nodes=32
#SBATCH --ntasks-per-node=48
#SBATCH --time=0:02:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3700


runFolder=$PWD"/MPI32_"$(date +"D_%Y_%m_%d_T_%I_%M_%S")
echo "Creating running folder" $runFolder

binary="/home/dcs/csrcnj/phd/OP-PIC/fempic_mpi/bin/mpi_hdf5"
echo "Using Binary" $binary

module load GCC/10.3.0  OpenMPI/4.1.1
module load PETSc/3.15.1
module load HDF5/1.12.1

# export I_MPI_PMI_LIBRARY=/usr/lib/x86_64-linux-gnu/libpmi.so
# export I_MPI_PMI_LIBRARY=/usr/lib64/pmix/lib/libpmi.so

# export OMP_PLACES=cores
# export OMP_PROC_BIND=close
# export OMP_NUM_THREADS=1

num_nodes=$SLURM_JOB_NUM_NODES

configFile="box_avon.param"
file=$PWD/$configFile

for run in 1 2 3; do
    for config in 1536000; do
    
        echo $file $config

        # ****************************************
        echo "Running MPI"

        folder=$runFolder/$config"_mpi"

        mkdir -p $folder
        cp $file $folder
        currentfilename=$folder/$configFile

        sed -i "s/STRING hdf_filename = <path_to_hdf5_mesh_file>/STRING hdf_filename = \/home\/dcs\/csrcnj\/phd\/box_mesh_gen\/hdf5_8\/box_${config}.hdf5/" ${currentfilename}

        srun $binary $currentfilename  | tee $folder/log_N${num_nodes}_C${config}_D10_R${run}.log;

        sed -i "s/REAL plasma_den     = 1e18/REAL plasma_den     = 1.3e18/" ${currentfilename}

        srun $binary $currentfilename  | tee $folder/log_N${num_nodes}_C${config}_D13_R${run}.log;

        # sed -i "s/REAL plasma_den     = 0.5e18/REAL plasma_den     = 1e17/" ${currentfilename}

        # srun $binary $currentfilename  | tee $folder/log_N${num_nodes}_C${config}_D1_R${run}.log;
    done
done

echo "simulation done"
