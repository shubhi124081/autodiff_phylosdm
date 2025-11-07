# ============================================================
# Temporal train/test indices per cluster & interval (R)
# ============================================================

rm(list = ls())

# ---------- User params ----------
root <- "~/phylo-sdms2"
dpath <- file.path(root, "raw_data/hummingbirdSA")
MAIN_DOMAIN <- "e2024" # FULL dataset file prefix
TEMPORAL_DOMAIN <- "n2005" # one of: n2005, n2012, n2015, n2018, FULL
NREPS <- 5
SEED_BASE <- 20251006

# CSV must have: species, cluster, interval_label, interval_start, interval_end, n_points
TEMPORAL_CSV <- file.path(dpath, "ebird_temporal_data_points_per_species_cluster_interval.csv")

# ---------- Helpers ----------
sample_n_safe <- function(vec, n, seed = NULL) {
    vec <- unique(as.integer(vec))
    if (length(vec) == 0L || n <= 0L) {
        return(integer(0))
    }
    if (!is.null(seed)) set.seed(seed)
    if (length(vec) <= n) {
        return(vec)
    }
    sample(vec, n)
}

# Build resolver name -> index from store$species_index_map (robust to layout)
make_species_resolver <- function(store) {
    if (!("species_index_map" %in% names(store))) {
        return(NULL)
    }
    m <- store$species_index_map
    # try to produce a named vector: names = species_name, values = integer index
    if (!is.null(names(m))) {
        # case A: names are indices (numeric-like), values are names
        if (all(suppressWarnings(!is.na(as.integer(names(m)))))) {
            idx <- as.integer(names(m))
            nm <- as.character(unname(m))
            return(setNames(idx, nm))
        } else {
            # case B: names are names, values are indices (or numeric-like)
            val <- suppressWarnings(as.integer(unname(m)))
            return(setNames(val, names(m)))
        }
    } else {
        # case C: vector of names, indices implied by position
        nm <- as.character(m)
        idx <- seq_along(nm)
        return(setNames(idx, nm))
    }
}

resolve_species_index <- function(sp_name, resolver) {
    if (!is.null(resolver) && sp_name %in% names(resolver)) {
        return(as.integer(resolver[[sp_name]]))
    }
    # allow numeric-like species ids in CSV
    maybe_idx <- suppressWarnings(as.integer(sp_name))
    if (!is.na(maybe_idx)) {
        return(maybe_idx)
    }
    stop(
        "Cannot resolve species '", sp_name, "' to an integer index. ",
        "Provide store$species_index_map or use numeric indices in the CSV."
    )
}

# ---------- Inputs ----------
# clusters
load(file.path(dpath, "spList.Rdata")) # loads spList (list: cluster -> character species names)
cluster_index <- seq_along(spList)

# temporal counts
temporal_n <- read.csv(TEMPORAL_CSV, stringsAsFactors = FALSE, check.names = FALSE)
req_cols <- c("species", "cluster", "interval_label", "interval_start", "interval_end", "n_points")
miss <- setdiff(req_cols, names(temporal_n))
if (length(miss)) stop("Temporal CSV missing columns: ", paste(miss, collapse = ", "))

# ---------- Main loop ----------
for (i in cluster_index) {
    cluster_name <- names(spList)[i]
    message("Cluster: ", cluster_name)

    # Load FULL store for this cluster
    store_path <- file.path(dpath, cluster_name, paste0(MAIN_DOMAIN, "_ALL_run_files.Rdata"))
    if (!file.exists(store_path)) {
        warning("Missing FULL store for cluster ", cluster_name, " at ", store_path, "; skipping.")
        next
    }
    load(store_path) # loads 'store'
    if (!exists("store")) {
        warning("No object 'store' in ", store_path, "; skipping.")
        next
    }

    # Basic store checks
    if (!is.list(store) || !all(c("y", "x") %in% names(store))) {
        warning("store missing 'y' or 'x' for cluster ", cluster_name, "; skipping.")
        next
    }
    y_full <- store$y

    need_y <- c("site", "species", "count")

    if (any(!need_y %in% names(y_full))) {
        warning("store$y missing cols; skipping ", cluster_name)
        next
    }

    # types
    y_full$site <- as.integer(y_full$site)
    y_full$species <- as.integer(y_full$species)
    y_full$count <- as.numeric(y_full$count)

    # CSV subset for this cluster and interval
    cfg <- subset(temporal_n, cluster == cluster_name & interval_label == TEMPORAL_DOMAIN)
    if (!nrow(cfg)) {
        warning("No rows in CSV for cluster ", cluster_name, " & interval ", TEMPORAL_DOMAIN, "; skipping.")
        next
    }

    # species name -> index resolver
    resolver <- make_species_resolver(store)

    # Build indices for NREPS
    list_of_indices <- vector("list", NREPS)

    for (rep_i in seq_len(NREPS)) {
        # per species storage
        train_by_sp <- list()
        test_by_sp <- list()

        # deterministic but rep-specific seeds
        rep_seed <- SEED_BASE + rep_i * 100000

        # iterate config rows (species)
        for (r in seq_len(nrow(cfg))) {
            sp_name <- as.character(cfg$species[r])
            sp_idx <- resolve_species_index(sp_name, resolver)
            n_take <- as.integer(cfg$n_points[r])
            yr_min <- as.integer(cfg$interval_start[r])
            yr_max <- as.integer(cfg$interval_end[r])

            sp_rows <- y_full[y_full$species == sp_idx, , drop = FALSE]
            if (!nrow(sp_rows)) {
                train_by_sp[[sp_name]] <- integer(0)
                test_by_sp[[sp_name]] <- integer(0)
                next
            }

            # presence / absence within interval (unique site ids)
            pres_sites <- unique(sp_rows$site[sp_rows$count > 0])
            abs_sites <- unique(sp_rows$site[sp_rows$count == 0])

            # sample independently
            pres_seed <- rep_seed + sp_idx * 19L
            abs_seed <- rep_seed + sp_idx * 17L
            pres_take <- sample_n_safe(pres_sites, n_take, seed = pres_seed)
            abs_take <- sample_n_safe(abs_sites, n_take, seed = abs_seed)

            train_sites <- sort(unique(c(pres_take, abs_take)))

            # test is "everything else for that species" (all years)
            all_sites_sp <- unique(sp_rows$site)
            test_sites <- sort(setdiff(all_sites_sp, train_sites))

            train_by_sp[[sp_name]] <- train_sites
            test_by_sp[[sp_name]] <- test_sites
        }

        # Flatten combined (unique)
        training_sites_combined <- sort(unique(unlist(train_by_sp, use.names = FALSE)))
        testing_sites_combined <- sort(unique(unlist(test_by_sp, use.names = FALSE)))

        list_of_indices[[rep_i]] <- list(
            training_sites = training_sites_combined,
            testing_sites = testing_sites_combined,
            training_by_species = train_by_sp,
            testing_by_species = test_by_sp
        )
    }

    names(list_of_indices) <- paste0("rep_", seq_len(NREPS))
    save_path <- file.path(dpath, cluster_name, paste0(TEMPORAL_DOMAIN, "_ALL_indices.Rdata"))
    save(list_of_indices, file = save_path)
    save(store, file = file.path(dpath, cluster_name, paste0(TEMPORAL_DOMAIN, "_ALL_run_files.Rdata")))
    message("  Saved: ", save_path)
}
