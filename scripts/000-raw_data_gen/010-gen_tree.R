# Set up --
rm(list = ls())

# Some libs
library(ape)
library(phytools)

# Filepaths & scripts
script.directory <- file.path(getwd(), "scripts")

# Read files in --

# Configuration parameters for data generation
config <- yaml::read_yaml(file.path(script.directory, "config.yaml"),
  as.named.list = TRUE
)
data_sim <- config$data_wrangling

# Commit config file now --
config$phylo$mean_root_vector <- yamlConvert(vec = TRUE)
config$phylo$covar_traits <- yamlConvert(vec = FALSE)
phylo <- config$phylo
tree <- fixTree(tree_time = phylo$tree_time)
tree <- checkTree()

print(paste0("# Tips : ", length(tree$tip.label)))
plot(tree)
tree$tip.label <- LETTERS[seq_len(length(tree$tip.label))]
# Save tree as a pdf
pdf(file = file.path(getwd(), "trees", paste0(phylo$tree_name, "Tree", ".pdf")))
plot(tree)
dev.off()

save(tree, file = file.path(
  getwd(), "trees",
  paste0(phylo$tree_name, ".Rdata")
))


########### For number of tips ############
if (!is.null(config$phylo$ntips)) {
  ntips <- config$phylo$ntips
  tree <- ape::rtree(n)
}

# Optional save

# save(tree, file = file.path(
#   getwd(), "trees",
#   paste0(phylo$tree_name, ".Rdata")
# ))
