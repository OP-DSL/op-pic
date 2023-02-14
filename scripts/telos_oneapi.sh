
. ~/.bashrc

module purge;
module load intel-oneapi-2021.3/intel-classic-2021.3;
module load intel-oneapi-2021.3/impi-classic-2021.3;
module load cuda/toolkit-10.2.89;
export LD_LIBRARY_PATH=/ext-home/zl/lib_install/petsc-3.18.2-oneapi-2021.3-release/lib:$LD_LIBRARY_PATH

# export LD_LIBRARY_PATH=/ext-home/zl/phd/OP-PIC/app_fempic_oppic/lib:/ext-home/zl/lib_install/petsc-3.18.2-oneapi-2021.3/lib:$LD_LIBRARY_PATH


