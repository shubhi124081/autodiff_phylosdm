library(data.table)

if (!interactive()) {
    args <- commandArgs(trailingOnly = TRUE)
    start <- args[1]
    stop <- args[2]
} else {
    start <- 1
    stop <- 1000
}
# master <- read.csv("~/phylo-sdms/phyloproj/raw_data/hummingbird_taxa_resolution.csv")
master <- read.csv("~/project/phylo-sdms/raw_data/hummingbird_taxa_resolution.csv")
ebird_list <- master$Master_sps
ebird_list <- ebird_list[!is.na(ebird_list)]
ebird_list <- gsub("_", " ", ebird_list)
setwd("/gpfs/gibbs/pi/jetz/data/species_datasets/occurrence/ebird/ebird_Apr2024/global_allcols/raw/all_files")
fls <- list.files()
ebrd_store <- list()
# length(fls)
for (iii in start:stop) {
    print(paste0(iii, " of ", length(fls)))
    ebrd <- read.csv(fls[iii], sep = ";", quote = "\"", as.is = TRUE)
    if (nrow(ebrd) == 0) {
        print("No data before")
        next
    }
    ebrd <- data.table(ebrd)
    ebrd <- ebrd[ebrd$taxonrank == "species", ]
    ebrd <- ebrd[ebrd$scientificname %in% ebird_list, ]
    ebrd <- ebrd[ebrd$approved == 1, ]
    ebrd <- ebrd[ebrd$all_species_reported == 1, ]
    ##
    ebrd$group_identifier <- ifelse(ebrd$group_identifier == "", NA, ebrd$group_identifier)
    ebrd$sampling_event_identifier <- ifelse(!is.na(ebrd$group_identifier), ebrd$group_identifier, ebrd$sampling_event_identifier)
    ebrd <- ebrd[(ebrd$samplingprotocol == "Traveling" |
        ebrd$samplingprotocol == "Stationary") &
        ebrd$effort_distance_km <= 1 &
        ebrd$duration_minutes >= 30]
    ##
    ebrd$year <- year(ebrd$eventDate)
    ebrd$day <- yday(ebrd$eventDate)
    ebrd$month <- month(ebrd$eventDate)
    ebrd <- ebrd[ebrd$year >= 2003 & ebrd$year <= 2024, ]
    if (nrow(ebrd) == 0) {
        print("No data after")
        next
    }
    ##
    ebrd <- unique(ebrd[, c(
        "scientificname",
        "sampling_event_identifier",
        "latitude",
        "longitude",
        "year",
        "month",
        "day",
        "samplingprotocol",
        "effort_area_ha",
        "effort_distance_km",
        "duration_minutes",
        "number_observers"
    )])
    ebrd_store[[iii]] <- ebrd
}

data <- do.call(rbind, ebrd_store)
setwd("~/project/ebird_processing")
save(data, file = paste0("ebird_start_", start, ".RData"))
