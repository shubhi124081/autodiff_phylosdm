# Purpose of script: Build bash script that will run 002_main which pulls
# together the models, dataset and config params to actually run the experiment
# This script will specify HPC specs, which model to run and which dataset to use

# Notes: 002_buildjobs.R and 002_main.R should be run in conjuction with
# each other. This script will specify dataset, model and all hpc params.
# The output of this script is a bash file (that can be submitted
# automatically if needed) that will call 002_main.R, with trailing args
# including dataset name and model name that 002_main.R will read in and run

# Set up ----
HPC <- Sys.getenv("HPC")
if (HPC == "FALSE") {
  root <- "~/phylo-sdms2"
} else {
  root <- "/vast/palmer/pi/jetz/ss4224/phylo-sdms2"
}

hpc_root_filepath <- "/vast/palmer/pi/jetz/ss4224/phylo-sdms2"

# Filepaths & scripts
scripts_directory <- file.path(root, "scripts")
data_directory <- file.path(root, "data")
job_directory <- file.path(root, "jobs")
log_directory <- file.path(root, "log")

##################### Read Config  #####################
# Configuration parameters for data generation
config <- yaml::read_yaml(file.path(scripts_directory, "000-config.yaml"),
  as.named.list = TRUE
)

phylo <- config$phylo
data_list <- config$data
data_sim <- config$data_sim
hpc_specs <- config$hpc_specs
buildLocal <- data_list$buildLocal
exp_root <- data_list$exp_root # The name/folder of datasets
exp_id <- data_list$exp_id # The name of the dataset
exp_name <- paste0(exp_root, "_", exp_id)
shortname <- paste0(exp_name, "_", data_list$cluster) # Used for naming the folders

Rscript_prefix <- "Rscript "
# Generally "Rscript " will work, but in case you need to point to another
# install of Rscript, change this

models <- config$stan_specs$models
datasets <- dir(file.path(data_directory, shortname))

expfiles <- expand.grid(x = models, y = datasets)
colnames(expfiles) <- c("model", "dataset")
expfiles$dataset <- gsub(".Rdata", "", expfiles$dataset)

header <- "
module load miniconda
conda activate brms"
shebang <- "#!/bin/bash"
dsq <- hpc_specs$dsq

# Two ways to submit - job array or job on loop
# way one
if (dsq) {
  modules <- "module load miniconda; conda activate brms; "
  Rscript_prefix <- "Rscript "

  command <- paste0(
    Rscript_prefix,
    file.path(hpc_root_filepath, "scripts", "023-main.R")
  )
  # Output file
  script1 <- gsub("\\.R", "", "023-main.R")
  exp_name1 <- shortname
  filename <- paste0("job_array_", exp_name1, "_for_", script1, ".txt")
  output_file <- file.path(job_directory, filename)
  nloop <- nrow(expfiles)

  # Open a connection to the file
  file_conn <- file(output_file, open = "wt")

  for (o in seq_len(nloop)) {
    command_fin <- paste(
      modules, command, expfiles[o, "model"],
      exp_name1, expfiles[o, "dataset"], "\n"
    )
    # writeLines(command_fin, output_file)
    cat(sprintf(command_fin), file = file_conn)
  }

  # Close the file connection
  close(file_conn)

  print(sprintf(
    "To submit, load dsq and create a job file with name %s",
    filename
  ))

  system(sprintf("job_name=%s", filename))
} else {
  if (buildLocal) {
    # If built local but run on the HPC, then write the hpc filepath into the
    # job script
    command <- paste0(
      Rscript_prefix,
      file.path(hpc_root_filepath, "scripts", "023-main.R"),
      " "
    )
  } else {
    command <- paste0(
      Rscript_prefix,
      file.path(scripts_directory, "023-main.R"),
      " "
    )
  }

  prefix <- "#SBATCH --"

  for (o in seq_len(nrow(expfiles))) {
    # jobname <- paste0(expfiles[o, "dataset"], "+", expfiles[o, "model"], ".sh")
    jobname_main <- paste0(shortname, "_", o, ".sh")
    job_file <- file(file.path(job_directory, jobname_main), "w")

    mid <- paste0(
      prefix, "job-name=", jobname_main, "\n",
      prefix, "partition=", hpc_specs$PARTITION, "\n",
      prefix, "out=", "log/slurm-%j.out", "\n",
      prefix, "time=", hpc_specs$TIME, "\n",
      prefix, "nodes=", hpc_specs$NODES, "\n",
      prefix, "ntasks=", hpc_specs$NTASKS, "\n",
      prefix, "cpus-per-task=", hpc_specs$CPUSPERTASK, "\n",
      prefix, "mem=", hpc_specs$MEM, "\n",
      prefix, "mail-user=", hpc_specs$MAIL_USER, "\n",
      prefix, "mail-type=", hpc_specs$MAIL_TYPE
    )

    command_fin <- paste(
      command, expfiles[o, "model"],
      shortname, expfiles[o, "dataset"]
    )

    writeLines(
      c(shebang, mid, "\n", header, command_fin),
      file.path(job_directory, jobname_main)
    )
    # writeLines(
    #   c(shebang, "\n", command_fin),
    #   file.path(job_directory, jobname_main)
    # )

    if (hpc_specs$AUTO_SUBMIT && user != "shubhi") {
      system(paste("sbatch", file.path(job_directory, jobname_main)))
    }

    if (!hpc_specs$AUTO_SUBMIT) {
      print(paste0(
        "To submit, run sbatch ", file.path(job_directory, jobname_main)
      ))
    }
  }
}
# write.csv(expfiles, file = paste0(log_directory, "/", shortname, ".csv"),
#           row.names = TRUE)
