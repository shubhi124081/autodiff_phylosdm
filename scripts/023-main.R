# Pull together stan models and generated data from bash script to run the
# experiment
# Set up ----
HPC <- Sys.getenv("HPC")
if (HPC == "FALSE") {
  root <- "~/phylo-sdms2"
} else {
  root <- "/vast/palmer/pi/jetz/ss4224/phylo-sdms2"
}

data_directory <- file.path(root, "data")
scripts_directory <- file.path(root, "scripts")
res_directory <- file.path(root, "res")
job_directory <- file.path(root, "jobs")
log_directory <- file.path(root, "log")
args <- commandArgs(trailingOnly = TRUE)

print(args)
model <- args[1]
shortname <- args[2] # shortname is the name of the folder
dataset <- args[3]

# Source all the stan models
source(file.path(scripts_directory, "012-stan_models.R"))

# Load data
filename <- paste0(dataset, ".Rdata")
load(file.path(data_directory, shortname, filename))

standata <- everything$data
standata$offset <- standata$offset + 1e-10 # avoid log(0)
standata$offset <- log(standata$offset)
config <- everything$config
data_sim <- config$data_sim
stan_specs <- config$stan_specs
code <- as.character(mostCommonlyUsed[model])

# TODO: Make sure there are no NAs in the data

# Wait, remove site
standata$X <- standata$X[, !colnames(standata$X) %in% "site"]
standata$K <- standata$K - 1

# Stan model
options(mc.cores = 4)
fit <- rstan::stan(
  model_code = code,
  data = standata,
  iter = stan_specs$iter,
  thin = stan_specs$thin,
  warmup = stan_specs$warmup,
  chains = stan_specs$chains,
  cores = stan_specs$cores
)
# Create a data directory and filepath
res_filepath <- file.path(res_directory, shortname)
if (!dir.exists(res_filepath)) {
  dir.create(res_filepath)
}
all <- list(
  fit = fit,
  config = config
  # data = standata
)
shortname_model <- paste0(dataset, "_", model)
save(all, file = file.path(res_filepath, paste0(shortname_model, ".Rdata")))

# Archive
# control = list(
#   adapt_delta = stan_specs$delta,
#   stepsize = stan_specs$stepsize,
#   stepsize_jitter = stan_specs$stepsize_jitter
# ),
# Chuck "useless" columns - file sizes are seriously bloated
# fit <- as.matrix(fit)
# b_ind <- grep("beta", colnames(fit))
# st_ind <- grep("st_devs", colnames(fit))
# y_ind <- grep("y_devs", colnames(fit))
# r_ind <- grep("rho", colnames(fit))
# a_ind <- grep("alpha", colnames(fit))
# indices <- c(b_ind, st_ind, y_ind, a_ind, r_ind)
# fit <- fit[, indices]

# model_log_name <- data_sim$experiment_name
# allfiles <- dir(log_directory)
# txtfiles <- allfiles[grepl("txt", allfiles)]
# txtfiles <- gsub(".txt", "", txtfiles)
# if (!(model_log_name %in% txtfiles)) {
#   write.table(" ", file = file.path(log_directory,
#               paste0(model_log_name, ".txt" )))
# }

# my_blank_line <- (" ")
# myline <- paste( model_log_name, "=> ", "model = ", model,
#                  "iteration = ", stan_specs$iter,
#                  "warmup = ", stan_specs$warmup, "chains = " ,stan_specs$chains,
#                  #"thin =", stan_specs$thin,
#                  "with dataset =", dataset)
# write(my_blank_line,
#       file = file.path(log_directory, paste0(model_log_name, ".txt" )),
#       append = TRUE)
# write(myline,                                            # Write new line to file
#       file =  file.path(log_directory, paste0(model_log_name, ".txt" )),
#       append = TRUE)
