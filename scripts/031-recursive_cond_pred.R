#!/usr/bin/env Rscript
# ============================================================
# 031-recursive_cond_pred.R
# Batch runner for 031-cond_pred.R across many experiment configs
# - Checkpoints: skips if all cluster outputs already exist
# - Logging: writes one rolling log file
# - Parallel: optional via MC_CORES env var
#
# USAGE
#   # 1) Use the built-in grid (edit GRID_* below):
#   Rscript scripts/031-recursive_cond_pred.R
#
#   # 2) Or provide a CSV manifest with columns:
#   #    EXP_ROOT,EXP_ID,FSP,DEF_LEV,REPNO,CLUSTERS
#   #    (CLUSTERS can be NA / blank to mean "all clusters")
#   Rscript scripts/031-recursive_cond_pred.R path/to/manifest.csv
# ============================================================
rm(list = ls())
suppressPackageStartupMessages({
    library(parallel)
    library(tools)
})

# ---------- Setup paths ----------
HPC <- Sys.getenv("HPC")
if (nzchar(HPC)) {
    root <- "~/phylo-sdms2"
} else {
    root <- "~/phylo-sdms2"
}
scripts_directory <- file.path(root, "scripts")
data_directory <- file.path(root, "data")
raw_directory <- file.path(root, "raw_data")
res_directory <- file.path(root, "res")
analysis_directory <- file.path(root, "analysis")
cond_pred_directory <- file.path(analysis_directory, "happy", "cond_pred") # created ad hoc later per EXP_ROOT

# Ensure logs dir exists
logs_dir <- file.path(analysis_directory, "logs")
if (!dir.exists(logs_dir)) dir.create(logs_dir, recursive = TRUE)
log_file <- file.path(logs_dir, sprintf("recursive_condpred_%s.log", format(Sys.time(), "%Y%m%d-%H%M%S")))

log_msg <- function(...) {
    msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste0(..., collapse = ""))
    cat(msg, "\n")
    cat(msg, "\n", file = log_file, append = TRUE)
}

# ---------- Discover clusters helper (matches logic inside 031-cond_pred.R) ----------
discover_clusters <- function(EXP_ROOT, EXP_ID) {
    # res/ contains folders/files named like "<EXP_ROOT>_<EXP_ID>_<CLUSTER>..."
    all_entries <- dir(res_directory, full.names = FALSE)
    pat <- sprintf("^%s_%s_", gsub("\\.", "\\\\.", EXP_ROOT), gsub("\\.", "\\\\.", EXP_ID))
    focal <- all_entries[grepl(pat, all_entries) & !grepl("\\.tar\\.gz$", all_entries)]
    # extract <CLUSTER> from "<ROOT>_<ID>_<CLUSTER>..."
    sub(paste0("^", EXP_ROOT, "_", EXP_ID, "_"), "", focal)
}

# ---------- Expected output checker (per combo) ----------
all_outputs_exist <- function(EXP_ROOT, EXP_ID, FSP, DEF_LEV, REPNO, clusters) {
    out_dir <- file.path(analysis_directory, EXP_ROOT, "cond_pred")
    if (!dir.exists(out_dir)) {
        return(FALSE)
    }
    exp_name <- paste0(EXP_ROOT, "_", EXP_ID)
    have <- vapply(clusters, function(cl) {
        shortname <- paste0(exp_name, "_", cl)
        filename <- paste0(shortname, "_", FSP, "_", DEF_LEV, "_", REPNO, "_cond_pred.Rdata")
        file.exists(file.path(out_dir, filename))
    }, logical(1))
    all(have)
}

# ---------- Command builder ----------
build_cmd <- function(EXP_ROOT, EXP_ID, FSP, DEF_LEV, REPNO, CLUSTERS) {
    # CLUSTERS: if NA/"" -> pass "NULL" so 031-cond_pred.R uses all clusters
    cl_arg <- if (is.na(CLUSTERS) || CLUSTERS == "") "NULL" else CLUSTERS
    sprintf(
        "Rscript %s/031-cond_pred.R %s %s %s %s %s %s",
        scripts_directory,
        shQuote(EXP_ROOT),
        shQuote(EXP_ID),
        shQuote(FSP),
        shQuote(DEF_LEV),
        shQuote(REPNO),
        shQuote(cl_arg)
    )
}

# ---------- Load spList to cross-check clusters if needed ----------
spList_path <- file.path(raw_directory, "hummingbirdSA", "spList.Rdata")
if (file.exists(spList_path)) {
    load(spList_path) # loads spList
} else {
    spList <- NULL
}

# ---------- Parameter grid ----------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 1 && file.exists(args[1])) {
    manifest <- read.csv(args[1], stringsAsFactors = FALSE, check.names = FALSE)
    req <- c("EXP_ROOT", "EXP_ID", "FSP", "DEF_LEV", "REPNO", "CLUSTERS")
    miss <- setdiff(req, names(manifest))
    if (length(miss)) stop("Manifest missing columns: ", paste(miss, collapse = ", "))
    grid <- manifest
} else {
    # EDIT these defaults as needed:
    # GRID_EXP_ROOT <- "happy"
    # GRID_EXP_ID <- c("AgCu2550100") # add more IDs here
    # GRID_FSP <- c("Aglaeactis_cupripennis") # or a specific "GENUS_SPECIES"
    # GRID_DEF_LEV <- c("y25", "y50", "y100") # e.g., e2024, y25, y50, y75
    # GRID_REPNO <- 1:5
    # GRID_CLUSTERS <- "Aglaeactis" # NA means "all clusters"
    # grid <- expand.grid(
    #     EXP_ROOT = GRID_EXP_ROOT,
    #     EXP_ID = GRID_EXP_ID,
    #     FSP = GRID_FSP,
    #     DEF_LEV = GRID_DEF_LEV,
    #     REPNO = GRID_REPNO,
    #     CLUSTERS = GRID_CLUSTERS,
    #     stringsAsFactors = FALSE
    # )

    # ---------- EDIT DEFAULT GRID HERE ----------

    # Define per-experiment specs
    # specs <- list(
    #     list(EXP_ID = "woopsie", FSP = "Coeligena_wilsoni", DEF_LEV = c("y1", "y2", "y5", "y10", "y25", "y50", "y100")),
    #     list(EXP_ID = "CoVi125102550", FSP = "Coeligena_violifer", DEF_LEV = c("y1", "y2", "y5", "y10", "y25", "y50")),
    #     list(EXP_ID = "CoTo125102550100", FSP = "Coeligena_torquata", DEF_LEV = c("y1", "y2", "y5", "y10", "y25", "y50", "y100")),
    #     list(EXP_ID = "CoPr125102550", FSP = "Coeligena_prunellei", DEF_LEV = c("y1", "y2", "y5", "y10", "y25", "y50")),
    #     list(EXP_ID = "CoPh12510", FSP = "Coeligena_phalerata", DEF_LEV = c("y1", "y2", "y5", "y10")),
    #     list(EXP_ID = "CoOr12", FSP = "Coeligena_orina", DEF_LEV = c("y1", "y2")),
    #     list(EXP_ID = "CoLu125102550", FSP = "Coeligena_lutetiae", DEF_LEV = c("y1", "y2", "y5", "y10", "y25", "y50")),
    #     list(EXP_ID = "CoIr125102550", FSP = "Coeligena_iris", DEF_LEV = c("y1", "y2", "y5", "y10", "y25", "y50")),
    #     list(EXP_ID = "CoHe125102550", FSP = "Coeligena_helianthea", DEF_LEV = c("y1", "y2", "y5", "y10", "y25", "y50")),
    #     list(EXP_ID = "CoCo125102550100", FSP = "Coeligena_coeligena", DEF_LEV = c("y1", "y2", "y5", "y10", "y25", "y50", "y100")),
    #     list(EXP_ID = "CoBo125102550", FSP = "Coeligena_bonapartei", DEF_LEV = c("y1", "y2", "y5", "y10", "y25", "y50"))
    # )
    specs <- list(
        list(CLUSTERS = "Abeillia"),
        list(CLUSTERS = "Adelomyia"),
        list(CLUSTERS = "Aglaeactis"),
        list(CLUSTERS = "Amazilia"),
        list(CLUSTERS = "Androdon"),
        list(CLUSTERS = "Anthracothorax"),
        list(CLUSTERS = "Archilochus"),
        list(CLUSTERS = "Boissonneaua"),
        list(CLUSTERS = "Chlorostilbon"),
        list(CLUSTERS = "Coeligena"),
        list(CLUSTERS = "Colibri"),
        list(CLUSTERS = "Discosura"),
        list(CLUSTERS = "Eriocnemis"),
        list(CLUSTERS = "Eugenes"),
        list(CLUSTERS = "Heliangelus"),
        list(CLUSTERS = "Heliodoxa"),
        list(CLUSTERS = "Phaethornis"),
        list(CLUSTERS = "Polytmus"),
        list(CLUSTERS = "Adelomyia2"),
        list(CLUSTERS = "Amazilia2"),
        list(CLUSTERS = "Amazilia3"),
        list(CLUSTERS = "Archilochus2"),
        list(CLUSTERS = "Phaethornis2")
    )


    # Shared params
    DEF_LEV <- "e2024"
    EXP_ID <- "woopsie"
    FSP <- "ALL"
    EXP_ROOT <- "happy"
    # CLUSTERS <- "Coeligena" # or NA_character_ to auto-discover
    REPNO <- 1
    TEMPORAL <- FALSE
    TT_LABEL <- "" # ignored when TEMPORAL=FALSE
    ART <- FALSE

    # Build rows
    #     rows <- lapply(specs, function(s) {
    #         expand.grid(
    #             EXP_ROOT = EXP_ROOT,
    #             EXP_ID = s$EXP_ID,
    #             FSP = s$FSP,
    #             DEF_LEV = s$DEF_LEV,
    #             REPNO = REPNO,
    #             CLUSTERS = CLUSTERS,
    #             TEMPORAL = TEMPORAL,
    #             TEMPORAL_TEST_DEF_LEV = TT_LABEL,
    #             ART = ART,
    #             stringsAsFactors = FALSE
    #         )
    #     })

    rows <- lapply(specs, function(s) {
        expand.grid(
            EXP_ROOT = EXP_ROOT,
            EXP_ID = EXP_ID,
            FSP = FSP,
            DEF_LEV = DEF_LEV,
            REPNO = REPNO,
            CLUSTERS = s$CLUSTERS,
            TEMPORAL = TEMPORAL,
            TEMPORAL_TEST_DEF_LEV = TT_LABEL,
            ART = ART,
            stringsAsFactors = FALSE
        )
    })
    grid <- do.call(rbind, rows)
    #     grid <- do.call(rbind, rows)
}

# Ensure output base directories exist
unique_exp_roots <- unique(grid$EXP_ROOT)
for (er in unique_exp_roots) {
    dir.create(file.path(analysis_directory, er, "cond_pred"), recursive = TRUE, showWarnings = FALSE)
}

# ---------- Runner (with checkpointing) ----------
run_one <- function(row_i) {
    r <- grid[row_i, , drop = FALSE]
    EXP_ROOT <- r$EXP_ROOT
    EXP_ID <- r$EXP_ID
    FSP <- r$FSP
    DEF_LEV <- r$DEF_LEV
    REPNO <- r$REPNO
    CLUSTERS <- r$CLUSTERS

    # Determine which clusters this combo would touch
    if (!is.na(CLUSTERS) && CLUSTERS != "") {
        clusters <- CLUSTERS
    } else {
        clusters <- discover_clusters(EXP_ROOT, EXP_ID)
        if (length(clusters) == 0) {
            log_msg(sprintf("SKIP (no clusters discovered): %s_%s", EXP_ROOT, EXP_ID))
            return(invisible(FALSE))
        }
    }

    # Checkpoint: skip if ALL outputs exist
    if (all_outputs_exist(EXP_ROOT, EXP_ID, FSP, DEF_LEV, REPNO, clusters)) {
        log_msg(sprintf(
            "SKIP (already complete): %s | %s | %s | %s | rep=%s",
            EXP_ROOT, EXP_ID, FSP, DEF_LEV, REPNO
        ))
        return(invisible(TRUE))
    }

    # Build and run the command
    cmd <- build_cmd(EXP_ROOT, EXP_ID, FSP, DEF_LEV, REPNO, if (is.na(CLUSTERS)) "" else CLUSTERS)
    log_msg(paste0("RUN  : ", cmd))

    # Capture stdout/stderr for the log
    out <- tryCatch(
        {
            res <- system(cmd, intern = TRUE, ignore.stderr = FALSE)
            cat(paste(res, collapse = "\n"), "\n", file = log_file, append = TRUE)
            TRUE
        },
        error = function(e) {
            log_msg(paste0("ERROR: ", conditionMessage(e)))
            FALSE
        }
    )

    # Re-check if complete now (partial failures will keep it “incomplete”)
    if (out) {
        if (all_outputs_exist(EXP_ROOT, EXP_ID, FSP, DEF_LEV, REPNO, clusters)) {
            log_msg(sprintf(
                "DONE : %s | %s | %s | %s | rep=%s",
                EXP_ROOT, EXP_ID, FSP, DEF_LEV, REPNO
            ))
        } else {
            log_msg(sprintf(
                "WARN : Post-run incomplete outputs for %s_%s (rep %s). Check logs.",
                EXP_ROOT, EXP_ID, REPNO
            ))
        }
    }
    invisible(out)
}

# ---------- Parallel or serial ----------
mc <- as.integer(Sys.getenv("MC_CORES", unset = "1"))
mc <- if (is.na(mc) || mc < 1) 1 else mc
log_msg(sprintf("Starting batch: %d jobs | MC_CORES=%d", nrow(grid), mc))

if (mc > 1 && .Platform$OS.type != "windows") {
    res <- mclapply(seq_len(nrow(grid)), run_one, mc.cores = mc)
} else {
    res <- lapply(seq_len(nrow(grid)), run_one)
}

ok_n <- sum(unlist(res), na.rm = TRUE)
log_msg(sprintf("Batch complete: %d/%d runs executed (skips counted as success).", ok_n, nrow(grid)))
