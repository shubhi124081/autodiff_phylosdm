#!/bin/bash

# Usage: ./download_species_files.sh species_list.txt remote_directory local_directory username@hpc_address

# Input arguments
SPECIES_LIST=$1       # Text file with species names, one per line
REMOTE_DIR=$2         # Directory on the HPC to search
LOCAL_DIR=$3          # Local directory to copy files to
HPC_USER=$4           # Username and HPC address (e.g., user@hpc.example.com)

# Check if all arguments are provided
if [ "$#" -ne 4 ]; then
  echo "Usage: $0 species_list.txt remote_directory local_directory username@hpc_address"
  exit 1
fi

# Create local directory if it doesn't exist
mkdir -p "$LOCAL_DIR"

# Loop through each species name in the species list
while IFS= read -r species_name; do
  echo "Searching for files matching: $species_name"

  # Use SSH to search for the file in the remote directory
  matching_files=$(ssh "$HPC_USER" "find '$REMOTE_DIR' -type f -name '$species_name*'")

  # Check if any files were found
  if [ -z "$matching_files" ]; then
    echo "No files found for $species_name"
  else
    echo "Files found for $species_name:"
    echo "$matching_files"

    # Use SCP to download each matching file
    for file in $matching_files; do
      echo "Downloading $file..."
      scp "$HPC_USER:$file" "$LOCAL_DIR"
    done
  fi
done < "$SPECIES_LIST"

echo "All done!"
