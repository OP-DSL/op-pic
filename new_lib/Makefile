
SRC = src
INC = src
LIB = lib
OBJ = obj

CPP        = mpicxx
CPPFLAGS   = -std=c++20 -O0 -g -fopenmp -fPIC #-Wall 
MPICPP	   = mpicxx
MPIFLAGS   = -std=c++20 -Wall -O3 -fopenmp -fPIC
NVCCFLAGS  = -O3 -m64 -Xptxas -dlcm=ca -Xptxas=-v -use_fast_math -arch=compute_60 -code=sm_60

CUDA_INSTALL_PATH = /opt/cuda/10.2.89/toolkit
CUDA_INC   = -I$(CUDA_INSTALL_PATH)/include
CUDA_LIB   = -L$(CUDA_INSTALL_PATH)/lib64

PETSC_INC  = # -I/opt/homebrew/Cellar/petsc/3.18.1/include
PETSC_LIB  = # -L/opt/homebrew/Cellar/petsc/3.18.1/lib -lpetsc

OPP_INC  = -I/opt/homebrew/Cellar/googletest/1.13.0/include
OPP_LIB  = -L/opt/homebrew/Cellar/googletest/1.13.0/lib -lgtest

.PHONY: clean mklib

ALL_TARGETS = clean test

all: $(ALL_TARGETS)

mklib:
	@mkdir -p $(LIB) $(OBJ)

oppic_lib: mklib
		$(CPP) $(CPPFLAGS) -o $(OBJ)/opp_structures.o -c $(SRC)/opp_structures.cpp
		$(CPP) $(CPPFLAGS) -o $(OBJ)/opp_simulation.o -c $(SRC)/opp_simulation.cpp
		$(CPP) $(CPPFLAGS) -o $(OBJ)/opp_params.o -c $(SRC)/opp_params.cpp

test: clean oppic_lib
		$(CPP) $(CPPFLAGS) -o $(OBJ)/test_opp_structures.o -c test/test_opp_structures.cpp $(OPP_INC)
		$(CPP) $(CPPFLAGS) -o $(OBJ)/test_opp_params.o -c test/test_opp_params.cpp $(OPP_INC)
		$(CPP) $(CPPFLAGS) -o $(OBJ)/test_main.o -c test/test_main.cpp $(OPP_INC)
		$(CPP) -std=c++20 -o test_main $(OBJ)/test_main.o $(OBJ)/test_opp_structures.o $(OBJ)/test_opp_params.o $(OBJ)/opp_structures.o $(OBJ)/opp_simulation.o $(OBJ)/opp_params.o $(OPP_INC) $(OPP_LIB)


clean:
		rm -f *.o *.d *.a
		rm -f $(OBJ)/*.o
		rm -f $(OBJ)/*~
		rm -f $(LIB)/*.a
		rm -f test_main
