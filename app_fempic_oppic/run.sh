#! /bin/bash

rm -rf files;
mkdir files;

rm -rf fempic

export OMP_NUM_THREADS=4
export OMP_PROC_BIND=close

# module load gnu-9.3.0/openmpi-4.0.4
# export LD_LIBRARY_PATH=/ext-home/zl/phd/OP-PIC/app_fempic_oppic:$LD_LIBRARY_PATH

# ./fempic_seq 

# ./fempic_genseq 

# ./fempic_openmp 

./fempic_genseq_wopet 

./fempic_openmp_wopet