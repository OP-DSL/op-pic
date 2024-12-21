#!/bin/bash
echo "Start date and time: $(date +"%Y-%m-%d %H:%M:%S")"

if [ $# -eq 0 ]; then
    echo "Usage: $0 <use_segmented_reduction::integer>"
    exit 1
fi

use_seg_red=$1
echo "Using Segmented reductions =" $use_seg_red

runFolder=$PWD"/MPI_SR"${use_seg_red}"_"$(date +"D_%Y_%m_%d_T_%I_%M_%S")

binpath=/cosma/home/do018/dc-lant1/phd/OP-PIC/app_fempic_cg/bin
binary=$binpath"/hip_mpi"

configFile="coarse.param"
file=/cosma/home/do018/dc-lant1/phd/OP-PIC/scripts/batch/fempic/cosma/$configFile

echo "Creating running folder -> " $runFolder
echo "Using Binary -> " $binary
echo "Config file -> " $file
echo "********************************************************"
cd $binpath
gitbranch=$(git branch | sed -n -e 's/^\* \(.*\)/\1/p')
gitcommit=$(git log -n 1 $gitbranch)
echo "Git branch " $gitbranch
echo "Git commit " $gitcommit
cd -
echo "********************************************************"

export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=close

plasma_den=1e18

for config in 48000 96000; do
    for run in 1 2 3; do
        for gpus in 2 4 8; do

            (( actual_config=config*gpus ))
            echo "Running MPI BLOCK " $actual_config $run $gpus

            folder=$runFolder/$config"_mpi"

            mkdir -p $folder
            cp $file $folder
            currentfilename=$folder/$configFile
      
            # escaped_folder="${folder//\//\\/}"
            sed -i "s|STRING global_mesh  = <path_to_mesh_files>/mesh.dat|STRING global_mesh  = /cosma/home/do018/dc-lant1/phd/box_mesh_gen/${actual_config}/mesh.dat|" ${currentfilename}
            sed -i "s|STRING inlet_mesh   = <path_to_mesh_files>/inlet.dat|STRING inlet_mesh   = /cosma/home/do018/dc-lant1/phd/box_mesh_gen/${actual_config}/inlet.dat|" ${currentfilename}
            sed -i "s|STRING wall_mesh    = <path_to_mesh_files>/wall.dat|STRING wall_mesh    = /cosma/home/do018/dc-lant1/phd/box_mesh_gen/${actual_config}/wall.dat|" ${currentfilename}

            if [ "$use_seg_red" -eq 1 ]; then
                sed -i "s/BOOL use_reg_red = false/BOOL use_reg_red = true/" ${currentfilename}
            fi

            sed -i "s/REAL plasma_den           = 1.0e18/REAL plasma_den     = ${plasma_den}/" ${currentfilename}
                        
            mpirun -np $gpus $binary $currentfilename > $folder/log_G${gpus}_M${actual_config}_D${plasma_den}_R${run}.log;
        done
    done
done

echo "End date and time: $(date +"%Y-%m-%d %H:%M:%S")"
