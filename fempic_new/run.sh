
rm -rf files/*;

export OMP_NUM_THREADS=1
export OMP_PROC_BIND=close

make PETSC=0 all

echo "BUILD DONE"

/ext-home/zl/phd/OP-PIC/fempic_new/bin/cuda /ext-home/zl/phd/OP-PIC/fempic_tests/configs/coarse_1g.param > xc1.log

echo "CUDA DONE"

/ext-home/zl/phd/OP-PIC/fempic_new/bin/seq /ext-home/zl/phd/OP-PIC/fempic_tests/configs/coarse_1gg.param > zz1.log

echo "SEQ DONE"