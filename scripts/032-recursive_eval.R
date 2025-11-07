#!/usr/bin/env Rscript
# ============================================================
# 032-recursive_eval.R
# Batch runner for scripts/032-eval.R across many configs
# - Discovers clusters from res/
# - Checkpoints (skips completed)
# - Logs runs
# - Parallel via MC_CORES
#
# Manifest CSV (optional) columns:
#   EXP_ROOT,EXP_ID,FSP,DEF_LEV,REPNO,CLUSTERS,TEMPORAL,TEMPORAL_TEST_DEF_LEV,ART
#   - CLUSTERS: NA/blank => discover all clusters
#   - TEMPORAL: TRUE/FALSE
#   - ART     : TRUE/FALSE
# ============================================================
rm(list = ls())
suppressPackageStartupMessages({
    library(parallel)
    library(tools)
})

# ---------- Paths ----------
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

# Logs
logs_dir <- file.path(analysis_directory, "logs")
if (!dir.exists(logs_dir)) dir.create(logs_dir, recursive = TRUE)
log_file <- file.path(logs_dir, sprintf(
    "recursive_eval_%s.log",
    format(Sys.time(), "%Y%m%d-%H%M%S")
))
log_msg <- function(...) {
    msg <- paste0(
        "[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
        paste0(..., collapse = "")
    )
    cat(msg, "\n")
    cat(msg, "\n", file = log_file, append = TRUE)
}

# ---------- Helpers ----------
# res/ entries look like "<EXP_ROOT>_<EXP_ID>_<CLUSTER>..."
discover_clusters <- function(EXP_ROOT, EXP_ID) {
    entries <- dir(res_directory, full.names = FALSE)
    pat <- sprintf(
        "^%s_%s_", gsub("\\.", "\\\\.", EXP_ROOT),
        gsub("\\.", "\\\\.", EXP_ID)
    )
    focal <- entries[grepl(pat, entries) & !grepl("\\.tar\\.gz$", entries)]
    sub(paste0("^", EXP_ROOT, "_", EXP_ID, "_"), "", focal)
}

# Expected eval filename for a (single) cluster
eval_filename_for <- function(EXP_ROOT, EXP_ID, FSP, DEF_LEV, REPNO, cluster) {
    exp_name <- paste0(EXP_ROOT, "_", EXP_ID)
    shortname <- paste0(exp_name, "_", cluster)
    paste0(shortname, "_", FSP, "_", DEF_LEV, "_", REPNO, "_eval.Rdata")
}

# Checkpoint: do all outputs already exist for this combo?
all_eval_outputs_exist <- function(EXP_ROOT, EXP_ID, FSP, DEF_LEV, REPNO, clusters) {
    out_dir <- file.path(analysis_directory, EXP_ROOT, "eval")
    if (!dir.exists(out_dir)) {
        return(FALSE)
    }
    all(vapply(clusters, function(cl) {
        file.exists(file.path(
            out_dir,
            eval_filename_for(EXP_ROOT, EXP_ID, FSP, DEF_LEV, REPNO, cl)
        ))
    }, logical(1)))
}

result_exists_for <- function(EXP_ROOT, EXP_ID, FSP, DEF_LEV, REPNO, cluster, model_name) {
    exp_name <- paste0(EXP_ROOT, "_", EXP_ID)
    shortname <- paste0(exp_name, "_", cluster)
    rfile <- file.path(
        res_directory, shortname,
        paste0(shortname, "_", FSP, "_", DEF_LEV, "_", REPNO, "_", model_name, ".Rdata")
    )
    file.exists(rfile)
}

# Build the Rscript command for 032-eval.R
build_cmd <- function(EXP_ROOT, EXP_ID, FSP, DEF_LEV, REPNO,
                      CLUSTERS, TEMPORAL, TEMPORAL_TEST_DEF_LEV, ART) {
    cl_arg <- if (is.na(CLUSTERS) || CLUSTERS == "") "NULL" else CLUSTERS
    tt_arg <- if (is.na(TEMPORAL_TEST_DEF_LEV)) "" else TEMPORAL_TEST_DEF_LEV
    sprintf(
        "Rscript %s/032-eval.R %s %s %s %s %s %s %s %s %s",
        scripts_directory, # 1
        shQuote(EXP_ROOT), # 2
        shQuote(EXP_ID), # 3
        shQuote(FSP), # 4
        shQuote(DEF_LEV), # 5
        shQuote(REPNO), # 6
        shQuote(cl_arg), # 7
        shQuote(as.logical(TEMPORAL)), # 8
        shQuote(tt_arg), # 9
        shQuote(as.logical(ART)) # 10
    )
}


# ---------- Parameter grid ----------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 1 && file.exists(args[1])) {
    manifest <- read.csv(args[1], stringsAsFactors = FALSE, check.names = FALSE)
    req <- c(
        "EXP_ROOT", "EXP_ID", "FSP", "DEF_LEV", "REPNO", "CLUSTERS",
        "TEMPORAL", "TEMPORAL_TEST_DEF_LEV", "ART"
    )
    miss <- setdiff(req, names(manifest))
    if (length(miss)) stop("Manifest missing columns: ", paste(miss, collapse = ", "))
    grid <- manifest
} else {
    # ---------- EDIT DEFAULT GRID HERE ----------

    # Define per-experiment specs

    # Define per-experiment specs
    # specs <- list(
    #     list(EXP_ID = "CoWi125102550100", FSP = "Coeligena_wilsoni", DEF_LEV = c("y1", "y2", "y5", "y10", "y25", "y50", "y100")),
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
        # list(CLUSTERS = "Discosura"),
        list(CLUSTERS = "Eriocnemis"),
        list(CLUSTERS = "Eugenes"),
        list(CLUSTERS = "Heliangelus"),
        list(CLUSTERS = "Heliodoxa"),
        list(CLUSTERS = "Phaethornis"),
        list(CLUSTERS = "Polytmus"),
        list(CLUSTERS = "Adelomyia2"),
        list(CLUSTERS = "Amazilia2"),
        # list(CLUSTERS = "Amazilia3"),
        # list(CLUSTERS = "Archilochus2"),
        list(CLUSTERS = "Phaethornis2")
    )

    # Shared params
    EXP_ROOT <- "happy"
    FSP <- "ALL"
    DEF_LEV <- "e2024"
    EXP_ID <- "woopsie"
    # CLUSTERS <- "Coeligena" # or NA_character_ to auto-discover
    REPNO <- 1
    TEMPORAL <- FALSE
    TT_LABEL <- "" # ignored when TEMPORAL=FALSE
    ART <- FALSE

    # # Build rows
    # rows <- lapply(specs, function(s) {
    #     expand.grid(
    #         EXP_ROOT = EXP_ROOT,
    #         EXP_ID = s$EXP_ID,
    #         FSP = s$FSP,
    #         DEF_LEV = s$DEF_LEV,
    #         REPNO = REPNO,
    #         CLUSTERS = CLUSTERS,
    #         TEMPORAL = TEMPORAL,
    #         TEMPORAL_TEST_DEF_LEV = TT_LABEL,
    #         ART = ART,
    #         stringsAsFactors = FALSE
    #     )
    # })

    # # Build rows
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

    # (Optional) Normalize types exactly once
    grid$REPNO <- as.character(grid$REPNO)
    grid$CLUSTERS <- ifelse(is.na(grid$CLUSTERS) | grid$CLUSTERS == "", NA_character_, grid$CLUSTERS)
    grid$TEMPORAL <- (as.logical(grid$TEMPORAL))
    grid$ART <- (as.logical(grid$ART))
    grid$TEMPORAL_TEST_DEF_LEV[is.na(grid$TEMPORAL_TEST_DEF_LEV)] <- ""
}

# Ensure per-root eval directories exist
for (er in unique(grid$EXP_ROOT)) {
    dir.create(file.path(analysis_directory, er, "eval"), recursive = TRUE, showWarnings = FALSE)
}

# ---------- Runner ----------
run_one <- function(i) {
    r <- grid[i, , drop = FALSE]
    EXP_ROOT <- r$EXP_ROOT
    EXP_ID <- r$EXP_ID
    FSP <- r$FSP
    DEF_LEV <- r$DEF_LEV
    REPNO <- r$REPNO
    CLUSTERS <- r$CLUSTERS
    TEMPORAL <- as.logical(r$TEMPORAL)
    TEMPORAL_TEST_DEF_LEV <- r$TEMPORAL_TEST_DEF_LEV
    ART <- as.logical(r$ART)

    # Which clusters are involved?
    clusters <- if (!is.na(CLUSTERS) && CLUSTERS != "") {
        CLUSTERS
    } else {
        discover_clusters(EXP_ROOT, EXP_ID)
    }

    if (length(clusters) == 0) {
        log_msg(sprintf("SKIP(no clusters): %s_%s", EXP_ROOT, EXP_ID))
        return(invisible(TRUE)) # treat as handled
    }

    # # Checkpoint: skip if ALL eval outputs exist for these clusters
    # if (all_eval_outputs_exist(EXP_ROOT, EXP_ID, FSP, DEF_LEV, REPNO, clusters)) {
    #     log_msg(sprintf(
    #         "SKIP(done): %s | %s | %s | %s | rep=%s",
    #         EXP_ROOT, EXP_ID, FSP, DEF_LEV, REPNO
    #     ))
    #     return(invisible(TRUE))
    # }

    # Determine clusters
    clusters <- if (!is.na(CLUSTERS) && CLUSTERS != "") CLUSTERS else discover_clusters(EXP_ROOT, EXP_ID)
    if (length(clusters) == 0) {
        log_msg(sprintf("SKIP(no clusters): %s_%s", EXP_ROOT, EXP_ID))
        return(invisible(TRUE))
    }

    # Require at least the baseline tag to exist (what 032-eval.R expects by default)
    baseline_ok <- vapply(
        clusters, function(cl) {
            result_exists_for(EXP_ROOT, EXP_ID, FSP, DEF_LEV, REPNO, cl, "LGCP_background")
        },
        logical(1)
    )

    if (!any(baseline_ok)) {
        log_msg(sprintf(
            "SKIP(no results): %s | %s | %s | %s | rep=%s",
            EXP_ROOT, EXP_ID, FSP, DEF_LEV, REPNO
        ))
        return(invisible(TRUE))
    }

    # Build and run command
    cmd <- build_cmd(
        EXP_ROOT, EXP_ID, FSP, DEF_LEV, REPNO,
        CLUSTERS, TEMPORAL, TEMPORAL_TEST_DEF_LEV, ART
    )
    log_msg(paste0("RUN : ", cmd))

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

    # Post-run check
    if (out) {
        if (all_eval_outputs_exist(EXP_ROOT, EXP_ID, FSP, DEF_LEV, REPNO, clusters)) {
            log_msg(sprintf(
                "DONE: %s | %s | %s | %s | rep=%s",
                EXP_ROOT, EXP_ID, FSP, DEF_LEV, REPNO
            ))
        } else {
            log_msg(sprintf(
                "WARN: Incomplete outputs after run: %s_%s (rep %s).",
                EXP_ROOT, EXP_ID, REPNO
            ))
        }
    }
    invisible(out)
}

# ---------- Parallel or serial ----------
mc <- as.integer(Sys.getenv("MC_CORES", unset = "1"))
mc <- if (is.na(mc) || mc < 1) 1 else mc

log_msg(sprintf("Starting eval batch: %d jobs | MC_CORES=%d", nrow(grid), mc))

if (mc > 1 && .Platform$OS.type != "windows") {
    res <- mclapply(seq_len(nrow(grid)), run_one, mc.cores = mc)
} else {
    res <- lapply(seq_len(nrow(grid)), run_one)
}

ok_n <- sum(unlist(res), na.rm = TRUE)
log_msg(sprintf(
    "Eval batch complete: %d/%d runs executed (some skips counted as success).",
    ok_n, nrow(grid)
))
