#!/bin/bash

rm -rf files;
mkdir files;

export OMP_NUM_THREADS=4
export OMP_PROC_BIND=close

# module load gnu-9.3.0/openmpi-4.0.4

CMDLINE1="bin/genseq /ext-home/zl/phd/OP-PIC/app_simpic_oppic/configs/system.param"
echo $CMDLINE1
$CMDLINE1

CMDLINE2="bin/omp /ext-home/zl/phd/OP-PIC/app_simpic_oppic/configs/system.param"
echo $CMDLINE2
$CMDLINE2

CMDLINE3="bin/seq /ext-home/zl/phd/OP-PIC/app_simpic_oppic/configs/system.param"
echo $CMDLINE3
$CMDLINE3

echo "Completed..."

