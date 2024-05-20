# OP-PIC_Configs

# Clone
git clone https://github.com/OP-DSL/OP-PIC.git

git clone https://github.com/ZamanLantra/OP-PIC_Configs.git

# Load Modules For Lumi

module load LUMI/23.09

module load partition/G

module load cray-hdf5-parallel

module load SCOTCH/6.1.3-cpeCray-23.09

export LD_LIBRARY_PATH=/project/project_465000743/INS_Code/petsc-install/lib:$LD_LIBRARY_PATH

export PETSC_INSTALL_PATH=/project/project_465000743/INS_Code/petsc-install

# Compile Lib

cd OP-PIC/lib_oppic

git checkout advec_mpi

export OPPIC_PATH=<path_to_lib_oppic_folder_in_OP-PIC_repo> 
OR ideally 
export OPPIC_PATH=$PWD

make hip_mpi

# Compile App

cd ../fempic_mpi/

make hip_mpi_hdf5

# Configure Slurm

cd <OP-PIC_Configs_folder>/fempic

python3 extract_tar.py

./update_slurm_files.sh <path_to_OP-PIC/fempic_mpi/bin_folder> <path_to_OP-PIC_Configs/fempic/mesh_files_folder>

# Run Slurm

cd lumi

sbatch fempic_mpi_1.sh

sbatch fempic_mpi_2.sh

sbatch fempic_mpi_4.sh

sbatch fempic_mpi_8.sh

sbatch fempic_mpi_16.sh

# Push the results to the same repo

Create a new branch in OP-PIC_Configs <Lumi_run_Sequence>

Add/Commit the files in fempic/lumi/MPI_D_*

push the changes