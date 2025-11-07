#!/bin/bash

# Define remote and local directories
remote_root="/gpfs/gibbs/project/jetz/ss4224/palmer_scratch/spatial_pred"
local_root="/Volumes/LaCie/phylo-sdms2/analysis/spatial_pred"

# Ensure local destination directory exists
mkdir -p "$local_root"

# Copy from HPC to local
scp mccleary:"$remote_root"/*hard_clipped.tif "$local_root"/

# Get list of subdirectories from remote root
#subdirs=$(ssh mccleary "ls $remote_root")

# Loop through each subdirectory
# for subdir in $subdirs; do
    # Define remote and local subdirectory paths
    # remote_subdir="$remote_root/$subdir/spatial_pred"
    # local_subdir="$local_root/$subdir/spatial_pred"

    # Create the local subdirectory if it does not exist
    # mkdir -p "$local_subdir"

    # Use scp to copy .tif files from remote to local directory
# scp mccleary:"$remote_subdir/pred_rast.tiff" "$local_subdir/"
# done


# To use rsync instead of scp, use 
$ rsync -avh --progress mccleary:~/palmer_scratch/spatial_pred/*hard_clipped.tif /Volumes/LaCie/phylo-sdms2/analysis/spatial_pred/