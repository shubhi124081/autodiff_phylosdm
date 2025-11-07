# ---- Inputs assumed present in parent frame ----
# raw, sps, tree_copy, cood_temporal, widenby, area_vect, epath, env_files,
# env_crs, eff_files, TEMPORAL, TEMPORAL_DOMAIN, dpath, cluster

# 1) Subset the data and tree -----------------------------------------------
stopifnot(all(sps %in% colnames(raw)))
y1 <- raw[, sps, drop = FALSE]

# Keep matching tips (by label safer than indices)
if (!all(sps %in% tree_copy$tip.label)) {
    # try to coerce label format if needed (e.g., spaces to underscores)
    alt <- gsub(" ", "_", sps)
    if (all(alt %in% tree_copy$tip.label)) sps <- alt
}
tree <- ape::keep.tip(tree_copy, tip = sps)

# 2) Cluster extent ----------------------------------------------------------
pres_rw <- which(rowSums(y1) > 0)
pres_cood <- cood_temporal[pres_rw, , drop = FALSE]
if (nrow(pres_cood) == 0) next

ex <- getExtentDf(pres_cood, NAMES = c("lat", "lon"))
ex <- terra::ext(ex + (c(-1, 1, -1, 1) * widenby)) # widen

# 3) Crop coords -------------------------------------------------------------
yc_crop <- cropDFCood(cbind(cood_temporal, y1), ex)

# 4) Empty raster at 1km -----------------------------------------------------
er <- createSampleRast(ex, area_vect, epath, env_files[1])

message("rasterizing...")
sps <- colnames(y1) # keep in sync with y1
NAMES <- c("x", "y")
data <- vector("list", length(sps))

# 5) Rasterize per species ---------------------------------------------------
# Loop over each species to rasterize occurrence points
for (i in seq_along(sps)) {
    sp <- sps[i]
    message(sp)
    # Find indices where the species is present
    ind1 <- which(yc_crop[, sp] == 1)
    if (length(ind1) == 0L) {
        # No presence points, store empty data.frame
        data[[i]] <- data.frame()
        next
    }

    # Extract coordinates and species column for presence points
    c1 <- yc_crop[ind1, c("x", "y", sp), drop = FALSE]
    # Create a spatial vector object from coordinates
    c1_vect <- terra::vect(c1[, NAMES, drop = FALSE], geom = NAMES)
    # Set CRS to match the empty raster (if possible)
    try(
        {
            terra::crs(c1_vect) <- terra::crs(er)
        },
        silent = TRUE
    )

    # Rasterize the points, summing occurrences per cell
    er1 <- terra::rasterize(c1_vect, er, fun = sum)
    # Convert raster to data.frame with coordinates
    er1_df <- terra::as.data.frame(er1, xy = TRUE)
    if (!nrow(er1_df)) {
        # No rasterized data, store empty data.frame
        data[[i]] <- data.frame()
        next
    }
    # Rename columns and add species label
    colnames(er1_df) <- c("x", "y", "sum")
    er1_df$sp <- sp
    # Store result for this species
    data[[i]] <- er1_df
}

names(data) <- sps
missing_species <- sps[vapply(data, nrow, 1L) == 0L]
if (length(missing_species)) {
    data <- data[!(names(data) %in% missing_species)]
    sps <- setdiff(sps, missing_species)
}

# 6) Wide site x species counts ---------------------------------------------
# (data is a list of data.frames with columns x,y,sum,sp)
if (length(data) == 0L) next
data_df <- data.table::rbindlist(data, use.names = TRUE, fill = TRUE)
data_sps <- data.table::dcast(
    data_df,
    x + y ~ sp,
    value.var = "sum",
    fun.aggregate = sum,
    fill = 0
)

# 7) Cap huge species (parameterize) ----------------------------------------
present_rows <- colSums(data_sps[, -(1:2), with = FALSE] > 0)
limit_flag <- present_rows > 5000L
if (any(limit_flag)) {
    keep_limit <- 5000L # <-- single source of truth
    sp_over <- names(present_rows)[limit_flag]
    for (sp in sp_over) {
        sp_rows <- which(data_sps[[sp]] > 0)
        # remove sites where other spp present to de-bias thinning
        other <- setdiff(names(data_sps)[-(1:2)], sp)
        cooccur_rows <- sp_rows[rowSums(data_sps[sp_rows, ..other] > 0) > 0]
        if (length(cooccur_rows)) data_sps <- data_sps[-cooccur_rows, ]
        sp_rows <- which(data_sps[[sp]] > 0) # refresh
        if (length(sp_rows) > keep_limit) {
            set.seed(42)
            drop_rows <- setdiff(sp_rows, sample(sp_rows, keep_limit))
            data_sps <- data_sps[-drop_rows, ]
        }
    }
}

# 8) Break into y / cood -----------------------------------------------------
y <- data_sps[, ..sps]
if (is.null(dim(y)) || nrow(y) <= 1L || sum(as.matrix(y)) <= 1L) next

cood <- as.data.frame(data_sps[, .(lon = x, lat = y)]) # keep order: lon,lat
cood <- cood[, c("lat", "lon")] # lat,lon for downstream

# 9) Annotate covariates -----------------------------------------------------
x <- annotateCoods(epath, env_files, env_crs, COOD = cood)
xx <- annotateCoods(
    "~/phylo-sdms2/raw_data/hummingbirdSA/effort_rasters",
    eff_files, env_crs, cood
)
colnames(xx) <- c("distance", "duration", "num_observers")
x <- cbind(x, xx)

# 10) Drop incomplete rows (fix >0) -----------------------------------------
rmv <- which(!complete.cases(x))
if (length(rmv) > 0L) {
    x <- x[-rmv, , drop = FALSE]
    y <- y[-rmv, , drop = FALSE]
    cood <- cood[-rmv, , drop = FALSE]
}

# 11) Handle non-finite safely ----------------------------------------------
x <- as.matrix(x)
nonfinite <- which(!is.finite(x))
if (length(nonfinite)) {
    message(sprintf("Replaced %d non-finite x entries with small eps.", length(nonfinite)))
    x[nonfinite] <- 1e-6
}
x <- process_x_matrix(as.data.frame(x))

if (is.null(nrow(y))) next

# 12) Synchronize rownames ---------------------------------------------------
rownames(x) <- rownames(cood) <- rownames(y) <- seq_len(nrow(y))

# 13) Background scaffold ----------------------------------------------------
J <- length(sps)
y$bg <- 0L
y <- as.data.frame(y) # Make sure it's a data frame before loop below

x_bg <- matrix(NA_real_, ncol = ncol(x), nrow = 1)
cood_bg <- matrix(NA_real_,
    ncol = ncol(cood), nrow = 1,
    dimnames = list(NULL, c("lon", "lat"))
)
y_bg <- matrix(NA_real_,
    ncol = ncol(y), nrow = 1,
    dimnames = list(NULL, colnames(y))
)
error_species <- integer(0)

# 14) Background per species -------------------------------------------------
for (species_idx in seq_len(J)) {
    message(sprintf("BG for species %d/%d", species_idx, J))
    tryCatch(
        {
            pres_idx <- which(y[, species_idx] > 0)
            if (length(pres_idx) < 2L) {
                warning("Not enough presence points, skipping BG")
                next
            }
            pres_cood <- cood[pres_idx, , drop = FALSE]
            ex <- getExtentDf(pres_cood, NAMES = c("lat", "lon"))

            bg_list <- createBgPoints(
                COOD = pres_cood, AREA_VECT = area_vect, EX = ex,
                NROW_Y = nrow(pres_cood), NCOL_Y = 1L, MULTIPLE = 1L, BUFFER = 50000
            )

            x_bg1 <- process_x_matrix(bg_list$x1)
            if (is.null(nrow(x_bg1)) || nrow(x_bg1) == 0L) {
                warning("No background points, skipping BG")
                error_species <- c(error_species, species_idx)
                next
            }
            # neutralize effort for BG
            for (nm in c("distance", "duration", "num_observers")) {
                if (nm %in% colnames(x_bg1)) x_bg1[, nm] <- 1
            }
            y_bg2 <- as.data.frame(matrix(0,
                nrow = nrow(x_bg1), ncol = ncol(y),
                dimnames = list(NULL, colnames(y_bg))
            ))
            y_bg2$bg <- species_idx

            x_bg <- rbind(x_bg, as.matrix(x_bg1))
            y_bg <- rbind(y_bg, as.matrix(y_bg2))
            cood_bg <- rbind(cood_bg, as.matrix(bg_list$bg_cood))
        },
        error = function(e) {
            warning(sprintf("Error processing species index %d: %s", species_idx, e$message))
            error_species <<- c(error_species, species_idx)
        }
    )
}

if (length(error_species)) {
    message("Species indices with BG errors: ")
    print(sort(unique(error_species)))
}

# 15) Drop incomplete BG scaffold rows --------------------------------------
rmv <- which(!complete.cases(x_bg) | !complete.cases(y_bg) | !complete.cases(cood_bg))
if (length(rmv) > 0L) {
    x_bg <- x_bg[-rmv, , drop = FALSE]
    y_bg <- y_bg[-rmv, , drop = FALSE]
    cood_bg <- cood_bg[-rmv, , drop = FALSE]
}

stopifnot(nrow(y_bg) == nrow(x_bg), nrow(y_bg) == nrow(cood_bg))

# Append BG
y <- rbind(as.data.frame(y), as.data.frame(y_bg))
x <- rbind(as.data.frame(x), as.data.frame(x_bg))
cood <- rbind(as.data.frame(cood), as.data.frame(cood_bg))

# 16) Site key ---------------------------------------------------------------
x$site <- y$site <- cood$site <- seq_len(nrow(y))

# 16.5) Scale columns ---------------------------------------------------------
# # scale environmental covariates (exclude effort vars and intercept)
nenv_vars <- c("Intercept", "site")
env_vars <- setdiff(colnames(x), nenv_vars)
xenv <- x[, env_vars, drop = FALSE]
rex_means <- apply(xenv, 2, mean, na.rm = TRUE)
rex_sds <- apply(xenv, 2, sd, na.rm = TRUE)
scales_df <- data.frame(variable = colnames(xenv), mean = rex_means, sd = rex_sds)

# Scale env vars
for (i in seq_len(ncol(xenv))) {
    varname <- colnames(xenv)[i]
    mean_val <- rex_means[i]
    sd_val <- rex_sds[i]
    if (!is.finite(sd_val) || sd_val == 0) {
        sd_val <- 1
        message(sprintf("Variable %s has non-finite or zero sd; setting sd to 1 to avoid division by zero.", varname))
    }
    xenv[, i] <- (xenv[, i] - mean_val) / sd_val
}

# Replace in x
if (length(env_vars) > 0L) {
    common <- intersect(colnames(x), colnames(xenv))
    if (length(common) == 0L) {
        warning("No matching env var columns found between x and xenv; nothing replaced.")
    } else {
        message("Replacing env var columns: ", paste(common, collapse = ", "))
        x[, common] <- xenv[, common, drop = FALSE]
    }
} else {
    message("No env_vars to replace.")
}

# Write out scales
# Save scales_df to file
write.csv(scales_df, file = file.path(dpath, paste0(cluster, "_", TEMPORAL_DOMAIN, "_ALL_env_scales.csv")), row.names = FALSE)

# 17) Fast long format (vectorized) -----------------------------------------
species_names <- sps
y_sponly <- as.matrix(y[, species_names, drop = FALSE])

nz <- which(y_sponly != 0, arr.ind = TRUE)
y_count <- data.frame(
    site    = nz[, 1],
    species = nz[, 2],
    count   = y_sponly[nz]
)

# add BG rows (count=0 for indicated species)
bg_rows <- which(y$bg > 0)
if (length(bg_rows)) {
    y_count <- rbind(
        y_count,
        data.frame(site = bg_rows, species = y$bg[bg_rows], count = 0)
    )
}

# 18) Tree: drop missing species (by label) ---------------------------------
if (length(missing_species)) {
    # ensure same naming convention as tree labels
    miss <- missing_species
    if (!all(miss %in% tree$tip.label)) {
        miss_alt <- gsub(" ", "_", miss)
        miss <- miss_alt[miss_alt %in% tree$tip.label]
    }
    if (length(miss)) tree <- ape::drop.tip(tree, miss)
}

# 19) Species index map (simple) --------------------------------------------
species_index_map <- setNames(species_names, seq_along(species_names))

# 19.5) If EXPERT_RANGE_OFFSET flag is on, need to add expert_range_offset --------
if (EXPERT_RANGE_OFFSET) {
    offset_path <- "/Volumes/LaCie/phylo-sdms2/analysis/soft_clips"
    if (!dir.exists(offset_path)) {
        warning("Expert range offset path not found. Please ensure Volumes/LaCie is connected.")
    }
    source("~/phylo-sdms2/scripts/000-raw_data_gen/009-expert_range_offset.R")
}

# 20) Save -------------------------------------------------------------------
message("saving...")
if (isTRUE(TEMPORAL)) {
    TEMPORAL_DOMAIN <- sub("e", "t", TEMPORAL_DOMAIN)
}
filename_train <- paste0(TEMPORAL_DOMAIN, "_ALL_run_files.Rdata")
dir.create(file.path(dpath, cluster), showWarnings = FALSE, recursive = TRUE)

# Save extent
ex_df <- data.frame(
    xmin = ex[1],
    xmax = ex[2],
    ymin = ex[3],
    ymax = ex[4]
)

write.csv(ex_df,
    file = file.path(dpath, paste0(cluster, "_", TEMPORAL_DOMAIN, "_ALL_extent.csv")), row.names = FALSE
)
store <- list(
    y = y_count, x = x, cood = cood, tree = tree,
    species_index_map = species_index_map
)
save(store, file = file.path(dpath, cluster, filename_train))

# 21) Clean up ---------------------------------------------------------------
rm(
    y1, ex, yc_crop, er, data_df, data_sps, y, cood, x, xx, rmv,
    y_count, present_rows, nz, bg_rows, store
)
gc()
