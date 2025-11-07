# EDA of eBird data for hummingbirds in South America
source("~/phylo-sdms2/scripts/000-raw_data_gen/000-raw_data_functions.R") # nolint
source("~/phylo-sdms2/scripts/000-phyloGenie_functions.R")

# Set-up
root <- "~/phylo-sdms2"
dpath <- file.path(root, "raw_data/hummingbirdSA")
epath <- "~/env"
TEMPORAL_DOMAIN <- "e2024" # Options: e2024, e2010, e2005, e2008
TEMPORAL <- FALSE # If TRUE, use temporal splits script

# Load spList
load(file.path(dpath, "spList.Rdata")) # loads obj spList

# Loop over clusters and store count_by_species for each
clusters <- names(spList)
count_list <- list()

for (CLUSTER in clusters) {
    raw_file <- file.path(
        dpath, CLUSTER, "e2024_ALL_run_files.Rdata"
    )
    contents <- load(raw_file)
    species_map <- store$species_index_map
    y <- store$y
    y$count <- ifelse(y$count > 0, 1, 0) # make binary
    count_by_species <- aggregate(y$count, by = list(y$species), FUN = sum)
    count_by_species$species_name <- species_map[as.character(count_by_species$Group.1)]
    count_by_species$cluster <- CLUSTER
    count_list[[CLUSTER]] <- count_by_species
}

all_counts <- do.call(rbind, count_list)
print(all_counts)

library(ggplot2)

# Plot for each cluster
for (CLUSTER in clusters) {
    df <- subset(all_counts, cluster == CLUSTER)
    p <- ggplot(df, aes(x = species_name, y = x)) +
        geom_bar(stat = "identity", fill = "steelblue") +
        geom_hline(yintercept = 200, linetype = "dotted", color = "red") +
        labs(
            title = paste("Species Data Points in Cluster:", CLUSTER),
            x = "Species",
            y = "Number of Data Points"
        ) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    print(p)
}
