#!/bin/bash

# This script updates binPath and hdfOriginalFolder in the .sh files in lumi folder using the inputs passed

if [ "$#" -ne 2 ]; then
    echo "Error: Please provide exactly two inputs."
    echo "Usage: $0 <path_to_bin> <path_to_mesh_files>"
    exit 1
fi

path_to_bin=$1
path_to_mesh_files=$2

echo "path_to_bin: $path_to_bin"
echo "path_to_mesh_files: $path_to_mesh_files"

escaped_path_to_bin="${path_to_bin//\//\\/}"
escaped_path_to_mesh_files="${path_to_mesh_files//\//\\/}"

slurm_file_folder=$PWD'/lumi'

# Iterate over all files in the folder
for file in "$slurm_file_folder"/*; do
    if [ -f "$file" ]; then
        echo "Updating file: $file"
        
        sed -i "s/binPath=<path_to_bin>/binPath=${escaped_path_to_bin}/" ${file}
        sed -i "s/hdfOriginalFolder=<path_to_mesh_files>/hdfOriginalFolder=${escaped_path_to_mesh_files}/" ${file}
    fi
done
