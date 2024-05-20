import os
import tarfile
from pathlib import Path

def untar_files(source_folder, destination_folder):
    # Make sure the destination folder exists
    Path(destination_folder).mkdir(parents=True, exist_ok=True)

    # Iterate over all files in the source folder
    for file_name in os.listdir(source_folder):
        file_path = os.path.join(source_folder, file_name)

        # Check if the file is a tar.gz file
        if file_name.endswith('.tar.gz'):
            print(f"Extracting {file_name}...")
            
            # Open the tar.gz file
            with tarfile.open(file_path, 'r:gz') as tar:
                # Extract all contents to the destination folder
                tar.extractall(destination_folder)

            print(f"{file_name} extracted successfully.\n")

if __name__ == "__main__":
    # Use the current working directory for source and destination folders
    source_folder = os.path.join(os.getcwd(), "mesh_files_tar")
    destination_folder = os.path.join(os.getcwd(), "mesh_files")

    untar_files(source_folder, destination_folder)
