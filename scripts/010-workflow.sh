#!/bin/bash

# Input parent directories below 
main=~/phylo-sdms2
main_mccleary=/vast/palmer/pi/jetz/ss4224/phylo-sdms2

# Run the tree generation script MANUALLY - 
# this was done on purpose so tree could be checked before experiment
# Run the data generation script 
big_string=$(Rscript $main/scripts/011-gen_data.R)
# echo $big_string
    
# Store experiment root name 
exp_name=$(echo $big_string | awk '{print $NF}')
echo $exp_name

Rscript $main/scripts/021-buildjobs.R

# Transfer files to hpc 
# scp $main/scripts/config.yaml mccleary:phylo-sdms/scripts
scp $main/scripts/023-main.R mccleary:$main_mccleary/scripts 
scp $main/scripts/012-stan_models.R mccleary:$main_mccleary/scripts 
scp $main/scripts/024-runjobs.R mccleary:$main_mccleary/scripts

# Copy .sh files that start with $exp_name
if ls $main/jobs/$exp_name*.sh 1> /dev/null 2>&1; then
  scp $main/jobs/$exp_name*.sh mccleary:$main_mccleary/jobs
fi

# Copy .txt files that contain $exp_name
if ls $main/jobs/*$exp_name*.txt 1> /dev/null 2>&1; then
  scp $main/jobs/*$exp_name*.txt mccleary:$main_mccleary/jobs
fi

# Transfer relevant dataset over 
cd $main/data
tar -zcvf $exp_name.tar.gz $exp_name
scp $main/data/$exp_name.tar.gz mccleary:$main_mccleary/data 
# NOTE: Need to be connected on Yale network for this 

# Open connection to hpc 
ssh mccleary

########## RUN A NON-DSQ JOB ON THE HPC ###########
# Now need to redefine variables 
main_mccleary=/vast/palmer/pi/jetz/ss4224/phylo-sdms2
late=$(find $main_mccleary/data -printf '%T+ %p\n' | sort -r | head -1)
y=${late%.tar.gz}
exp_name=${y##*/}
echo $exp_name

# Unzip files 
cd $main_mccleary/data
tar -zxvf $main_mccleary/data/$exp_name.tar.gz

# Important to change back to root directory otherwise sbatch freaks out 
# # (because it can't write the logfiles)
cd ~ 

# # Set up conda environment 
module load miniconda 
conda activate brms

# Run buildjob files - if autosubmit, they should just submit 
# Rscript $main_mccleary/scripts/021-buildjobs.R
cd $main_mccleary
Rscript $main_mccleary/scripts/024-runjobs.R $exp_name

########## RUN A DSQ JOB ON THE HPC ###########

# Now need to redefine variables 
main_mccleary=/vast/palmer/pi/jetz/ss4224/phylo-sdms2
late=$(find $main_mccleary/data -printf '%T+ %p\n' | sort -r | head -1)
y=${late%.tar.gz}
exp_name=${y##*/}
echo $exp_name

# Unzip files 
cd $main_mccleary/data
tar -zxvf $main_mccleary/data/$exp_name.tar.gz

# Load modules 
module load miniconda 
conda activate brms
module load dSQ

# Get job name 
job_name=$(ls -t $main_mccleary/jobs/*.txt | head -n 1 | xargs -n 1 basename)
echo $job_name
cd ~
dsq --job-file $main_mccleary/jobs/$job_name --mem-per-cpu 100g -t 24:00:00 --mail-type ALL --mail-type ALL --partition day --mail-user shubhi.sharma@yale.edu --out=log/slurm-%j.out

# save job id, useful commands below! 

########## POST-PROCESSING ON HPC NOT DSQ ###########
nano ~/project/phylo-sdms/scripts/config_post_process.yaml 
main_mccleary=~/project/phylo-sdms
module load miniconda 
conda activate phylo-sdms_wf
cd $main_mccleary 
Rscript $main_mccleary/scripts/040-build_jobs_pp.R

######## THE NEW POST-PROCESSING WF ON HPC W DSQ FOR SPATIAL PRED ########## 
main=~/phylo-sdms/phyloproj
scp $main/scripts/config_post_process.yaml mccleary:project/phylo-sdms/scripts
# On HPC 
main_mccleary=~/project/phylo-sdms
module load miniconda
conda activate phylo-sdms_wf
export TMPDIR=~/palmer_scratch 
module load dSQ 
Rscript $main_mccleary/scripts/040-build_jobs_pp.R 
# manual step: create $job_name - this is the txt file name printed from previous step
dsq --job-file $main_mccleary/jobs/$job_name --mem-per-cpu 100g -t 24:00:00 --mail-type ALL --mail-type ALL --partition day --mail-user shubhi.sharma@yale.edu

# save job id 
jobid=37276702
# ids 37587726, 37590291

# useful commands 
dsqa -j $jobid
cat job_37276702_status.tsv

########## IF YOU WANT TO ZIP UP FILES ########## 

# Zip up experiment files (assuming variables still exist, otherwise redefine)
main_mccleary=/vast/palmer/pi/jetz/ss4224/phylo-sdms2
cd $main_mccleary/res
tar -zcvf $exp_name.tar.gz $exp_name

# Run this to grab all the result files 
scp mccleary:$main_mccleary/res/$exp_name.tar.gz $main/res

########## IF YOU WANT TO UNZIP UP FILES ########## 
# Unzip 
cd $main/res
tar -zxvf $main/res/$exp_name.tar.gz
# say "beep beep" 

######### POST-PROCESSING ###########
# Need to run 025-conditional_prediction 
# Rscript $main/scripts/030-post_process_main.R $exp_name 
# say "I've finished running, your code is amazing. Good job"

# NOTE: To bring down a range map from gibbs
# sp_name=X  
mkdir $main/expert_ranges/$sp_name 
scp mccleary:/gpfs/gibbs/pi/jetz/data/species_datasets/rangemaps/birds/jetz_maps/jetz_maps_v2020/shapes/$sp_name/$sp_name.shp $main/expert_ranges/$sp_name
scp mccleary:/gpfs/gibbs/pi/jetz/data/species_datasets/rangemaps/birds/jetz_maps/jetz_maps_v2020/shapes/$sp_name/$sp_name.dbf $main/expert_ranges/$sp_name
scp mccleary:/gpfs/gibbs/pi/jetz/data/species_datasets/rangemaps/birds/jetz_maps/jetz_maps_v2020/shapes/$sp_name/$sp_name.prj $main/expert_ranges/$sp_name
scp mccleary:/gpfs/gibbs/pi/jetz/data/species_datasets/rangemaps/birds/jetz_maps/jetz_maps_v2020/shapes/$sp_name/$sp_name.shx $main/expert_ranges/$sp_name

# NOTE: To delete files from an experiment after it is archived, use this 
ls | grep "grumpy" | xargs rm -rf

# Note note: this will match everything so if you don't want to delete for e.g.
# grumpy2, match grumpy_

# Note note note: this command takes a while 
# Count with 

ls -1 | wc -l

# Run an interactive session 
srun -p devel --mem 28g --pty /bin/bash

# Figure out memory usage after a job has been completed (need jobid) 
sacct -j 37498895 --format=JobID,JobName,MaxRSS,MaxVMSize


# CoDel = 56617116
# CoSer = 56617338
# CoCor = 56617641
# CoTha = 56617686


#### Automate zipping of files ####
main_mccleary=/vast/palmer/pi/jetz/ss4224/phylo-sdms2
cd $main_mccleary/res

# loop over all happy_n05_* directories
for exp_name in happy_woopsie_*; do
  if [ -d "$exp_name" ]; then
    echo "Zipping $exp_name ..."
    tar -zcvf "${exp_name}.tar.gz" "$exp_name"
  fi
done


#### Automate copying and unzipping of files ####
main=~/phylo-sdms2
main_mccleary=/vast/palmer/pi/jetz/ss4224/phylo-sdms2

cd $main/res

# copy all tarballs down
scp mccleary:$main_mccleary/res/happy_woopsie_*.tar.gz .

# unzip each one
for f in happy_n15_*.tar.gz; do
  echo "Unzipping $f ..."
  tar -zxvf "$f"
done

echo "beep beep 🚗"


#### Automate copying and unzipping of files ####
main=~/phylo-sdms2
main_mccleary=/vast/palmer/scratch/jetz/ss4224/spatial_pred

cd $main/analysis/happy/spatial_pred

# copy all tarballs down
scp mccleary:$main_mccleary/happy_BigBkg4_Aglaeactis_*.tif .

# unzip each one
#for f in happy_n15_*.tar.gz; do
#  echo "Unzipping $f ..."
#  tar -zxvf "$f"
#done

echo "beep beep 🚗"
