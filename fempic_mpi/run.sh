
rm -rf files/*;

export OMP_NUM_THREADS=1
export OMP_PROC_BIND=close
export OMP_PLACES=cores

make all

echo "BUILD DONE"

/ext-home/zl/phd/OP-PIC/fempic_new/bin/cuda /ext-home/zl/phd/OP-PIC/fempic_tests/configs/coarse_3.param

echo "CUDA DONE"

/ext-home/zl/phd/OP-PIC/fempic_new/bin/seq /ext-home/zl/phd/OP-PIC/fempic_tests/configs/coarse_3.param

echo "SEQ DONE"

export OMP_NUM_THREADS=48
/ext-home/zl/phd/OP-PIC/fempic_new/bin/omp /ext-home/zl/phd/OP-PIC/fempic_tests/configs/coarse_3.param

mpirun -np 48 valgrind --log-file=vg%p --track-origins=yes bin/mpi /ext-home/zl/phd/OP-PIC/scripts/fempic_tests/configs/coarse_5.param

echo "OMP DONE"

# bin/cuda /ext-home/zl/phd/OP-PIC/fempic_tests/configs/medium_1.param -vec_type cuda -mat_type aijcusparse
# https://www.youtube.com/watch?v=jiOWnG_Kk-U

# watch -n 0.1 /opt/rocm-5.4.3/bin/rocm-smi

# bin/cuda /ext-home/zl/phd/OP-PIC/scripts/fempic_tests/configs/coarse_5.param
# nvprof --analysis-metrics -o nv_dy_sh.nvvp bin/cuda /ext-home/zl/phd/OP-PIC/scripts/fempic_tests/configs/coarse_5.param


# module load forge/arm-forge-22.1.3
# export ALLINEA_FORCE_LICENCE_FILE=/opt/arm/forge/licences/licences_x86_64/License_Temp
# export ALLINEA_FORCE_LICENCE_FILE=/opt/arm/forge/licences/licences_x86_64/License_Temp

# mpirun -np 40 vtune -collect performance-snapshot –r vtune_mpi40_ps -- bin/mpi
# vtune-gui vtune_mpi8.telos/
# 	performance-snapshot
# 	hotspots
# 	anomaly-detection
# 	memory-consumption
# 	uarch-exploration
# 	memory-access
# 	threading
# 	hpc-performance
# 	io
# 	gpu-offload
# 	gpu-hotspots
# 	fpga-interaction
# 	system-overview
# 	graphics-rendering
# 	cpugpu-concurrency
# 	vpp
# 	tsx-exploration
# 	tsx-hotspots
# 	sgx-hotspots
# 	uarch-exploration-0

# advixe-cl -collect roofline -project-dir MyResults -- bin/omp /ext-home/zl/phd/OP-PIC/scripts/fempic_tests/configs/coarse_5.param

# /opt/cuda/10.2.89/toolkit/nsight-compute-2019.5.0/nv-nsight-cu-cli

/opt/cuda/10.2.89/toolkit/bin/nsys profile --stats=true --output=nv_compute2 bin/cuda /ext-home/zl/phd/OP-PIC/scripts/fempic_tests/configs/coarse_5.param

/opt/cuda/10.2.89/toolkit/nsight-compute-2019.5.0/nv-nsight-cu-cli --export nvcompute2 bin/cuda /ext-home/zl/phd/OP-PIC/scripts/fempic_tests/configs/coarse_5.param

# spack load intel-oneapi-compilers@2023.1.0
# spack load parmetis@4.0.3
# spack load hipsycl@0.9.4
# spack load hdf5@1.14.1-2
# spack load cmake@3.26.3
# spack load intel-mpi@2019.10.317
# spack load cube@4.8
# spack load scalasca@2.6.1

# export SCOREP_ENABLE_TRACING=true

# spack load cube@4.8
# spack load otf2@3.0
# spack load cubew@4.8
# spack load binutils@2.40

# export LD_LIBRARY_PATH=/ext-home/zl/lib_setup/spack/opt/spack/linux-debian10-skylake_avx512/gcc-8.3.0/otf2-3.0-g6ehsuayrfll2mx3wwnykkjesxmporav/lib:$LD_LIBRARY_PATH
# export LD_LIBRARY_PATH=/ext-home/zl/lib_setup/spack/opt/spack/linux-debian10-skylake_avx512/gcc-8.3.0/cube-4.8-nconzs6ywwpjwzgldwgfegtfpiewsyax/lib:$LD_LIBRARY_PATH
# export LD_LIBRARY_PATH=/ext-home/zl/lib_setup/spack/opt/spack/linux-debian10-cascadelake/intel-2021.3.0/binutils-2.40-lpfwg7ugeghftq7acfrc43pxxitfgttg/lib:$LD_LIBRARY_PATH
# export LD_LIBRARY_PATH=/ext-home/zl/lib_setup/spack/opt/spack/linux-debian10-skylake_avx512/gcc-8.3.0/binutils-2.40-phaf3hafohunayotjydvol7zixq6frgr/lib:$LD_LIBRARY_PATH

# scalasca -analyze mpirun -n 2 bin/mpi


# AVON ===============================================================

module load GCC/10.3.0  OpenMPI/4.1.1
module load PETSc/3.15.1
export PETSC_INSTALL_PATH=/scrtp/avon/eb/software/PETSc/3.15.1-foss-2021a
export OPPIC_PATH=/home/dcs/csrcnj/phd/OP-PIC/lib_oppic
export NVCCFLAGS_ADD='-gencode arch=compute_75,code=sm_75'
module load CUDA

# For MPI
# salloc -p hmem --nodes=1 --ntasks-per-node=48 --cpus-per-task=1 --mem-per-cpu=3700 --time=03:00:00 

# mpirun -np 8 bin/mpi /home/dcs/csrcnj/phd/OP-PIC/scripts/fempic_tests/configs/coarse_5.param

# salloc -p compute --nodes=2 --ntasks-per-node=48 --mem-per-cpu=3700 --time=03:00:00 --job-name=csrcnj
salloc -p gpu --nodes=1 --ntasks-per-node=48 --mem-per-cpu=3700 --gres=gpu:quadro_rtx_6000:3 --time=06:00:00 --job-name=csrcnj

mpirun -np 2 bin/cuda_mpi /home/dcs/csrcnj/phd/OP-PIC/scripts/fempic_tests/configs/coarse_2.param

mpirun -np 12 -hostfile hosts bin/cuda_mpi /home/dcs/csrcnj/phd/OP-PIC/scripts/fempic_tests/configs/coarse_2.param

# ===============================================================

# BEDE ===============================================================

module load gcc/12.2 openmpi/4.0.5 cuda/12.0.1 openblas/0.3.10
export PETSC_INSTALL_PATH=/users/csrcnl/lib_install/petsc-3.20.1_rel
export OPPIC_PATH=/users/csrcnl/phd/OP-PIC/lib_oppic
export NVCCFLAGS_ADD='-gencode arch=compute_70,code=sm_70'

# ===============================================================
