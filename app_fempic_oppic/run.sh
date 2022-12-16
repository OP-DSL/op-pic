#! /bin/bash

rm -rf files;
mkdir files;

rm -rf fempic

export OMP_NUM_THREADS=4
export OMP_PROC_BIND=close

# ./fempic_seq 

# ./fempic_genseq 

# ./fempic_openmp 

./fempic_genseq_wopet 

# ./fempic_openmp_wopet