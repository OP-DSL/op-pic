# OP-PIC FemPIC

In fempic.h, funcs_from_fempic.cpp and field_solver_xx.cpp, most of the original functions from https://www.particleincell.com/2015/fem-pic/ remains. These will be used mainly for loading the mesh.
field_solver_original.cpp contains the original code used in fempic which has Matrix computations (linear solver) written explicitly using c++.
field_solver_petsc.cpp contains all the other functions similar to field_solver_original.cpp, however the linear solve is done using sparse matrix KSP Solver of petsc.

The content in the omp and seq folders and any other file with suffix _op (fempic_op.cpp) should be auto generated once the source to source translator is developed.
For now, only the auto generated code works (both sequential and OpenMP)

Therefore, OP-PIC expects the user to write only fempic.cpp and kernels.h files.


