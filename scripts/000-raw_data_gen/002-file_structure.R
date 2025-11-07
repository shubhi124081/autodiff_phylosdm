# Define root directory
# Set up ----
HPC <- Sys.getenv("HPC")
root <- ifelse(HPC == FALSE, "~/phylo-sdms2", "~/project/phylo-sdms2")

# Define paths for analysis and raw data
analysis_dir <- file.path(root, "analysis")
raw_data_dir <- file.path(root, "raw_data/hummingbirdSA")

# Define experiment root directory
exp_root <- "happy"

# Load species list data
species_list_file <- file.path(raw_data_dir, "spList.RData")
load(species_list_file)

# Define other set of sub-directories
sub_dirs <- c("eval", "cond_pred", "figures", "spatial_pred")

# Create main exp_root directory
dir.create(file.path(analysis_dir, exp_root), recursive = TRUE)

for (species in names(spList)) {
    species_dir <- file.path(analysis_dir, exp_root, species)
    dir.create(species_dir)
    for (sub_dir in sub_dirs) {
        dir.create(file.path(species_dir, sub_dir))
    }
}
