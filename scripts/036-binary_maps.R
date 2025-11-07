# This script uses pre-computed thresholds to create binary maps from spatial prediction rasters.
# It applies expert range masks where available.
# It loops over clusters, model tags, and species to generate the binary maps.
# =============================== SETUP =======================================

rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
    library(terra)
    library(pROC)
})

# ----------------------------- Configuration ------------------------------

HPC <- Sys.getenv("HPC")
DATA_SAVER <- TRUE
EXP_ROOT <- "happy"

# Roots and output locations (local vs HPC)
if (HPC == "FALSE") {
    root <- "~/phylo-sdms2"
    message("Running locally")
    tempdir <- function() "/Volumes/LaCie/temp"
    output_dir <- file.path("/Volumes/LaCie/phylo-sdms2", "analysis", EXP_ROOT, "spatial_pred")
} else {
    root <- "/vast/palmer/pi/jetz/ss4224/phylo-sdms2"
    message("Running on HPC")
    tempdir <- function() "/vast/palmer/scratch/jetz/ss4224"
    output_dir <- "/vast/palmer/scratch/jetz/ss4224/spatial_pred"
}

# Prefer external tempdir if present (local)
if (dir.exists("/Volumes/LaCie/temp")) {
    terra::terraOptions(tempdir = "/Volumes/LaCie/temp")
} else {
    message("External temp drive not found. Using default tempdir.")
}

message("tempdir(): ", tempdir())
message("output_dir: ", output_dir)

# Handy cleaner (optional)
cleanup_terra_temp <- function(path = terra::terraOptions()$tempdir, days_old = 0) {
    files <- list.files(path, full.names = TRUE)
    if (!length(files)) {
        return(invisible())
    }
    old_files <- files[file.info(files)$mtime < (Sys.time() - days_old * 86400)]
    if (length(old_files)) unlink(old_files, recursive = TRUE)
    message(length(old_files), " files removed from terra tempdir.")
}

# ----------------------------- Paths --------------------------------------

scripts_directory <- file.path(root, "scripts")
data_directory <- file.path(root, "data")
raw_directory <- file.path(root, "raw_data")
res_directory <- file.path(root, "res")
analysis_directory <- file.path(root, "analysis")
expert_directory <- file.path(root, "expert_ranges")
eval_directory <- "~/phylo-sdms2/analysis/happy/eval"


pred_rast_directory <- if (DATA_SAVER) {
    if (HPC == "FALSE" || HPC == "") {
        file.path("/Volumes/LaCie/phylo-sdms2", "analysis", "pred_rast")
    } else {
        file.path(root, "analysis", "pred_rast")
    }
} else {
    file.path(analysis_directory, "pred_rast")
}

# ------------------------ Experiment parameters ---------------------------

if (interactive()) {
    EXP_ROOT <- "happy" # Root name for experiment
    EXP_ID <- "woopsie" # Experiment ID
    FSP <- "ALL"
    DEF_LEV <- "e2024"
    REPNO <- "1"
    CLUSTER <- "Coeligena" # one cluster name or NULL for all
    TEMPORAL <- FALSE
    ART <- FALSE
} else {
    args <- commandArgs(trailingOnly = TRUE)
    EXP_ROOT <- args[1]
    EXP_ID <- args[2]
    FSP <- args[3]
    DEF_LEV <- args[4]
    REPNO <- args[5]
    CLUSTER <- args[6]
    TEMPORAL <- as.logical(args[7])
    ART <- as.logical(args[8])
}
DEBUG <- FALSE
if (DEBUG) {
    print("DEBUG mode is ON")
    output_dir <- "~/Downloads/flash"
}

# Derived identifiers (per-cluster later)
exp_name <- paste0(EXP_ROOT, "_", EXP_ID)

# ------------------------------ Helpers -----------------------------------

# Decide model label by TAG (string → filename part)
model_name_for_tag <- function(tag) {
    if (tag == "NOPHYLO_BASE") "LGCP_nophylo_background" else "LGCP_background"
}

# Prediction raster path helper
pred_file_path <- function(out_dir, shortname, species, def_lev, repno, tag, trailing = "relprob") {
    file.path(out_dir, paste0(shortname, "_", species, "_", def_lev, "_", repno, "_", tag, "_", trailing, ".tif"))
}
# ------------------------------ Load base data -----------------------------

source(file.path(scripts_directory, "000-phyloGenie_functions.R"))
source(file.path(scripts_directory, "032-eval_functions.R"))

# species list by cluster
load(file.path(raw_directory, "hummingbirdSA", "spList.Rdata")) # loads spList

cluster_list <- if (CLUSTER == "NULL" || is.null(CLUSTER)) names(spList) else CLUSTER

TAGS <- c("PHYLO_BASE", "PHYLO_ALL", "NOPHYLO_BASE")

# =============================== MAIN LOOPS ===================================
for (cluster in cluster_list) {
    fspALL <- spList[[cluster]]
    J <- length(fspALL)
    shortname <- paste0(exp_name, "_", cluster)
    data_name <- paste0(shortname, "_", FSP, "_", DEF_LEV, "_", REPNO)

    for (tag in TAGS) {
        print(tag)
        species_indices <- if (ART) {
            which(fspALL == FSP)
        } else {
            seq_len(J)
        }

        # Read in thresholds
        for (species_idx in species_indices) {
            fsp <- fspALL[species_idx]

            tryCatch(
                {
                    # ------------------------ Threshold & write binaries ------------------------

                    pred_in <- file.path(
                        output_dir,
                        paste0(shortname, "_", fsp, "_", DEF_LEV, "_", REPNO, "_", tag, "_", "relprob", ".tif")
                    )
                    if (!file.exists(pred_in)) {
                        warning("Binary write skipped; missing: ", pred_in)
                        next
                    }

                    spat <- terra::rast(pred_in)

                    # Get threshold file
                    thr_file <- read.csv(file.path(
                        analysis_directory, EXP_ROOT, "eval", "thresholds",
                        paste0(shortname, "_", fsp, "_", DEF_LEV, "_", REPNO, "_", tag, "_thresholds.csv")
                    ))

                    # Wait, let's try loading the eval files to get thresholds and metrics
                    filename2 <- paste0(shortname, "_", FSP, "_", DEF_LEV, "_", REPNO, "_eval.Rdata")
                    contents <- load(file.path(eval_directory, filename2)) # loads obj results_df
                    results_df_sp <- results_df[results_df$species == fsp, ]
                    results_df_sp_tag <- results_df_sp[results_df_sp$Tag == tag, ]

                    # thresholding
                    thr <- thr_file$threshold[species_idx]
                    # thr <- results_df_sp_tag$threshold
                    if (is.na(thr)) thr <- 0
                    bin <- terra::ifel(spat > thr, 1, 0)

                    # hard clip to expert range
                    range_shp <- file.path(expert_directory, fsp, paste0(fsp, ".shp"))
                    if (file.exists(range_shp)) {
                        range <- terra::vect(range_shp)
                        # Buffer by 500km (it's in meters)
                        # range_buffer <- terra::buffer(range, width = 50000)
                        ex <- terra::ext(range)
                        ex <- terra::ext(ex[1] - 4, ex[2] + 4, ex[3] - 4, ex[4] + 4)
                        bin <- terra::crop(bin, ex)
                        terra::plot(bin, main = paste("Binary Map -", fsp, "-", tag, "-", DEF_LEV))
                        terra::lines(range)
                    } else {
                        warning("Expert range missing; writing unmasked binary for ", fsp)
                    }

                    out_bin <- file.path(
                        output_dir,
                        paste0(shortname, "_", fsp, "_", DEF_LEV, "_", REPNO, "_", tag, "_BINARY_THRESHOLD.tif")
                    )
                    terra::writeRaster(bin, filename = out_bin, overwrite = TRUE)
                },
                error = function(e) {
                    warning("Error processing species '", fsp, "' (index ", species_idx, "): ", conditionMessage(e))
                    # continue to next species
                }
            )
        } # species loop
    } # TAG loop
} # cluster loop
