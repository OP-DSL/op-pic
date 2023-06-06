cmake_minimum_required(VERSION 3.12)
project(OPPIC)

# 1) intel-2021.4/intel-oneapi-2021.4   2) intel-2021.4/impi-oneapi-2021.4
# set(CMAKE_C_COMPILER "/opt/intel/oneapi/2021.3/mpi/latest/bin/mpiicc")
# set(CMAKE_CXX_COMPILER "/opt/intel/oneapi/2021.3/mpi/latest/bin/mpiicpc")

# set(CMAKE_VERBOSE_MAKEFILE ON)

# set(PETSC_DIR /ext-home/zl/lib_install/petsc-3.18.3-oneapi-2021.3-mpi_telos)
set(CMAKE_CXX_FLAGS "-std=c++11 -fPIC -Wwrite-strings")

set(USE_PETSC ${PETSC} CACHE STRING "Description of the argument")

if(USE_PETSC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_PETSC")
endif()
# set(OpenMP_CXX_FLAGS "-fopenmp")
# set(OpenMP_CXX_LIB_NAMES "omp")

find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)
# find_package(ParMETIS REQUIRED)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

set(CORE_SOURCES
    src/core/oppic_lib_core.cpp 
    src/core/oppic_util.cpp 
    src/core/trace.cpp 
    src/core/opp_params.cpp 
)

# Add library target for sequential version
add_library(opp_seq STATIC 
    ${CORE_SOURCES}
    src/seq/oppic_seq.cpp
)

# Add library target for OpenMP version
add_library(opp_omp STATIC 
    ${CORE_SOURCES}
    src/omp/oppic_omp.cpp
)
target_compile_options(opp_omp PRIVATE ${OpenMP_CXX_FLAGS})

# Add library target for MPI version
add_library(opp_mpi STATIC 
    ${CORE_SOURCES}
    src/mpi/opp_mpi.cpp
    src/mpi/opp_mpi_halo.cpp
    src/mpi/opp_mpi_partition.cpp
    src/mpi/opp_mpi_utils.cpp
    src/mpi/opp_mpi_particle_comm.cpp 
    src/mpi/opp_mpi_double_ind_reducs.cpp
)
target_compile_options(opp_mpi PRIVATE 
    ${MPI_CXX_COMPILE_OPTIONS}
)

# # Add library target for CUDA version
# find_package(CUDA REQUIRED)
# add_library(lib_oppic_cuda STATIC oppic.cu)
# target_include_directories(lib_oppic_cuda PRIVATE ${CUDA_INCLUDE_DIRS})
# target_link_libraries(lib_oppic_cuda PRIVATE ${CUDA_LIBRARIES})

# # Add executable target
# add_executable(oppic main.cpp)

# # Link libraries to executable
# target_link_libraries(oppic PRIVATE lib_oppic_seq)  # Sequential version
# target_link_libraries(oppic PRIVATE lib_oppic_omp)  # OpenMP version
# target_link_libraries(oppic PRIVATE lib_oppic_mpi)  # MPI version
# target_link_libraries(oppic PRIVATE lib_oppic_cuda)  # CUDA version
