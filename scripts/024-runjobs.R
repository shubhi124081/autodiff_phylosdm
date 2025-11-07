# This is a stupid script that I wrote because I need to submit multiple jobs
# at once and don't want to invest time in learning how to
# do it in bash
data_directory <- file.path(getwd(), "data")
scripts_directory <- file.path(getwd(), "scripts")
res_directory <- file.path(getwd(), "res")
job_directory <- file.path(getwd(), "jobs")
log_directory <- file.path(getwd(), "log")

args <- commandArgs(trailingOnly = TRUE)
dataset_main <- args[1]
all_jobs <- dir(job_directory)
my_jobs <- all_jobs[grepl(paste0(dataset_main, "_"), all_jobs)]
n_jobs <- length(my_jobs)

for (i in 1:n_jobs) {
  system(paste("sbatch", file.path(job_directory, my_jobs[i])))
}
