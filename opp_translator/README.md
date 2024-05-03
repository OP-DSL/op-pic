# OPP Code Generator

### Requirements

- python >= 3.8
- libclang

source the required file in OP-PIC/source/ directory

change directory to OP-PIC/opp_translator

run setup_venv.sh -- one time process

build the lib file (for now, till it gets automated) using the make file in OP-PIC/opp_lib folder 
[e.g., make seq; make mpi; make cuda; make cuda_mpi; etc...]

python3 $OPP_INSTALL_PATH/../opp_translator/translator -v -I $OPP_INSTALL_PATH/../opp_lib/include/ --file_paths cabana.cpp