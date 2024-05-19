# OP-PIC SimPIC

In simpic.h and funcs_from_simpic.cpp, most of the original functions from https://bitbucket.org/lecad-peg/simpic/src/master/ remains. These will be used mainly for loading the mesh and particles.

The algorithm contains the original code used in simpic which has Matrix computations written explicitly using c++, which cannot be easily converted to kernel calls.
