
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

echo "OMP DONE"

# bin/cuda /ext-home/zl/phd/OP-PIC/fempic_tests/configs/medium_1.param -vec_type cuda -mat_type aijcusparse
# https://www.youtube.com/watch?v=jiOWnG_Kk-U

# watch -n 0.1 /opt/rocm-5.4.3/bin/rocm-smi


# module load forge/arm-forge-22.1.3
# export ALLINEA_FORCE_LICENCE_FILE=/opt/arm/forge/licences/licences_x86_64/License_Temp
# export ALLINEA_FORCE_LICENCE_FILE=/opt/arm/forge/licences/licences_x86_64/License_Temp

