rm(list = ls())

# Example clade Colibri
root <- "~/phylo-sdms2"
raw_directory <- file.path(root, "raw_data", "hummingbirdSA")

# Load the run files
CLUSTER <- "Coeligena"
raw_file <- file.path(
    raw_directory, CLUSTER, "e2024_ALL_run_files.Rdata"
)
contents <- load(raw_file)
tree <- store$tree
species_map <- store$species_index_map
y <- store$y
y$count <- ifelse(y$count > 0, 1, 0) # make binary
count_by_species <- aggregate(y$count, by = list(y$species), FUN = sum) # check presences per species
count_by_species$species_name <- species_map[as.character(count_by_species$Group.1)]
print(count_by_species)
x <- store$x
cood <- store$cood

# EDA plot
# spno <- 2
# y_species <- y[y$species == spno, ]
# cood_species <- cood[cood$site %in% y_species$site, ]
# contents <- load("~/phylo-sdms2/raw_data/world.Rdata") # world
# yc <- cbind(cood_species, y_species)
# yc <- yc[, c("lon", "lat", "count")]
# colnames(yc) <- c("x", "y", "sp")
# makeOccurrenceMap(yc, "sp")

# Basically we'll create a whole new set of indices
sps <- tree$tip.label
DEF_LEV <- "y100"
DEF <- 100
NREP <- 5

for (j in seq_len(length(sps))) {
    fsps <- sps[j]
    cat("Processing species:", fsps, "\n")

    # Create indices for fsps
    fsps_index <- as.numeric(names(species_map)[which(species_map == fsps)])

    # Need to make each species data-deficient
    y_species <- y[y$species == fsps_index, ]
    y_others <- y[y$species != fsps_index, ]

    # Sample sites
    n <- nrow(y_species)
    n_pres <- length(which(y_species$count > 0))
    if ((n_pres * 0.75) < DEF) {
        cat("  Skipping species:", fsps, "because it has only", n_pres, "presences\n")
        next
    }
    list_of_indices <- list()

    for (i in seq_len(NREP)) {
        cat("  Replicate:", i, "of", NREP, "\n")
        pres_sites <- y_species[y_species$count > 0, ]
        abs_sites <- y_species[y_species$count == 0, ]

        # Sample 5 presences and 5 absences
        train_pres_idx <- sample(seq_len(nrow(pres_sites)), size = DEF, replace = FALSE)
        train_abs_idx <- sample(seq_len(nrow(abs_sites)), size = DEF, replace = FALSE)

        train_pres_sites <- pres_sites[train_pres_idx, ]
        train_abs_sites <- abs_sites[train_abs_idx, ]

        # Combine training sites (5 presences + 5 absences) and all y_others
        # Clean y_others
        y_others <- y_others[!y_others$site %in% train_pres_sites$site, ]
        # Take focal presence sites for the focal species out
        y_tmp <- y[y$site %in% y_others$site, ]
        y_tmp_species <- y_tmp[y_tmp$species == fsps_index, ]
        if (nrow(y_tmp_species) > 0) {
            rmv_sites <- y_tmp_species$site
            y_others <- y_others[!y_others$site %in% rmv_sites, ]
        }

        train_all <- rbind(train_pres_sites, train_abs_sites, y_others)

        # Test sites: the rest of y_species not in training
        test_pres_sites <- pres_sites[-train_pres_idx, , drop = FALSE]
        test_abs_sites <- abs_sites[-train_abs_idx, , drop = FALSE]
        test_sites <- rbind(test_pres_sites, test_abs_sites)

        list_of_indices[[i]] <- list(
            training_sites = train_all$site,
            testing_sites = test_sites$site
        )
    }

    names(list_of_indices) <- paste0("rep_", seq_len(NREP))
    out_file <- file.path(raw_directory, CLUSTER, paste0(DEF_LEV, "_", fsps, "_indices.Rdata"))
    save(list_of_indices, file = out_file)
    cat("  Saved indices to:", out_file, "\n")
    run_file_out <- file.path(
        raw_directory, CLUSTER, paste0(DEF_LEV, "_", fsps, "_run_files.Rdata")
    )
    y <- rbind(train_all, test_sites)
    x <- x[x$site %in% y$site, , drop = FALSE]
    cood <- cood[cood$site %in% y$site, , drop = FALSE]
    store <- list(
        y = y,
        x = x,
        cood = cood,
        species_index_map = species_map,
        tree = tree
    )
    save(store, file = run_file_out)
    cat("  Saved run_files to:", run_file_out, "\n")
}

# NOTE: training_sites and testing_sites are duplicated across species which is
# not ideal but we will just need to make sure "1" means "1". Background points
# are generated so they're not an issue

rm(list = ls())
contents <- load("~/phylo-sdms2/raw_data/hummingbirdSA/Aglaeactis/y1_Aglaeactis_castelnaudii_indices.Rdata") # list_of_indices
contents <- load("~/phylo-sdms2/raw_data/hummingbirdSA/Aglaeactis/y1_Aglaeactis_castelnaudii_run_files.Rdata") # store

y <- store$y
training_idx <- list_of_indices$rep_1$training_sites

y_train <- y[y$site %in% training_idx, ]
y_train_species <- y_train[y_train$species == 3, ] # store$species_index_map

cat("Number of training points for Aglaeactis castelnaudii:", nrow(y_train_species), "\n")
cat("Number of presences:", sum(y_train_species$count > 0), "\n")
