# Clear workspace
rm(list = ls())

# Libraries and functions (only what we use here)
source("~/phylo-sdms2/scripts/000-raw_data_gen/000-raw_data_functions.R") # nolint
source("~/phylo-sdms2/scripts/000-phyloGenie_functions.R")

# Set-up
root <- "~/phylo-sdms2"
dpath <- file.path(root, "raw_data/hummingbirdSA")

# Load datasets
load(file.path(dpath, "ebird_raw-2024.Rdata")) # loads obj: raw
colnames(raw)[colnames(raw) %in% c("latitude", "longitude")] <- c("lat", "lon")
colnames(raw) <- gsub(" ", "_", colnames(raw))

# Keep only rows with essential fields (avoid dropping due to other NA columns)
req_cols <- c("lat", "lon", "year")
missing_req <- setdiff(req_cols, colnames(raw))
if (length(missing_req)) stop("Raw data missing required columns: ", paste(missing_req, collapse = ", "))
raw <- raw[complete.cases(raw[, req_cols, drop = FALSE]), ]

# Load clusters
load(file.path(dpath, "spList.Rdata")) # loads obj: spList

# Year intervals
year_intervals <- list(
    n2005 = c(2003, 2005),
    n2012 = c(2003, 2012),
    n2015 = c(2003, 2015),
    n2018 = c(2003, 2018),
    FULL  = c(2003, 2024)
)

# Results accumulator
res <- data.frame(
    species = character(),
    cluster = character(),
    interval_label = character(),
    interval_start = integer(),
    interval_end = integer(),
    n_points = integer(),
    stringsAsFactors = FALSE
)

# Optional: list to store per-(cluster,species,interval) quick summaries
# store <- list()

# Loop over clusters
for (i in seq_along(spList)) {
    cluster_name <- names(spList)[i]
    fspALL <- spList[[i]] # character vector of species col names for this cluster

    # loop species in this cluster
    for (j in seq_along(fspALL)) {
        fsps <- as.character(fspALL[j])

        # guard: skip if species column doesn't exist
        if (!fsps %in% colnames(raw)) {
            warning(sprintf("Species column '%s' not found in raw; skipping.", fsps))
            next
        }

        # base species subset (species indicator + year)
        raw_sub0 <- raw[, c(fsps, "year"), drop = FALSE]

        # loop year intervals (independent filters each time)
        for (k in seq_along(year_intervals)) {
            label <- names(year_intervals)[k]
            yrs <- year_intervals[[k]]
            yr_min <- yrs[1]
            yr_max <- yrs[2]

            raw_int <- raw_sub0[raw_sub0$year >= yr_min & raw_sub0$year <= yr_max, , drop = FALSE]

            # number of points = sum of 1s in that species column within interval
            n <- if (nrow(raw_int)) sum(raw_int[[fsps]], na.rm = TRUE) else 0L

            # append a row
            res <- rbind(
                res,
                data.frame(
                    species = fsps,
                    cluster = cluster_name,
                    interval_label = label,
                    interval_start = as.integer(yr_min),
                    interval_end = as.integer(yr_max),
                    n_points = as.integer(n),
                    stringsAsFactors = FALSE
                )
            )

            # optional store entry keyed by cluster/species/interval label
            # key <- paste(cluster_name, fsps, label, sep = "__")
            # store[[key]] <- list(n_points = n, n_rows = nrow(raw_int))
        }
    }

    # memory valve (if raw is gigantic and you truly need to free)
    # rm(raw_sub0); gc()
}

# Save results
write.csv(res, file.path(dpath, "ebird_temporal_data_points_per_species_cluster_interval.csv"), row.names = FALSE)
