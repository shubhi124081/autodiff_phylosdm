# This script is going to be used to create the new ebird dataset
# The final ebird datast is compiled on the hpc using 000-compile_ebird_hpc.R using MOL ebird data
# The format of that data is record - species name - lat - long - date - time - protocol - area - distance - duration - observers
# This needs to be reformatted into record - lat - long etc - species names with 0 for non-detection and 1 for detection
# Set up ----
rm(list = ls())
t1 <- Sys.time()
HPC <- Sys.getenv("HPC")
if (HPC == FALSE) {
    contents <- load("~/phylo-sdms/phyloproj/raw_data/hummingbirdSA/ebird_final-2024.Rdata") # loads object data_all
} else {
    contents <- load("~/project/ebird_processing/ebird_final-2024.Rdata") # loads object data_all
}

if (!interactive()) {
    args <- commandArgs(trailingOnly = TRUE)
    start <- args[1]
    stop <- args[2]
    num <- args[3]
} else {
    start <- 1800001
    stop <- 2395468
    num <- 10
}

# Load compiled on the hpc
data_all$occurrence <- rep(1, nrow(data_all))

# Start and stop
all_brbs <- unique(data_all$scientificname)

# Dummy dataset
data_all2 <- data_all[start:stop, ]

# Takes too long with the whole dataset
# library(reshape2)

# # Create a wide format of the data with species as columns
# data_wide <- dcast(data_all2, sampling_event_identifier + latitude + longitude + year + month + day + samplingprotocol +
#     effort_area_ha + effort_distance_km + duration_minutes + duration_minutes ~ scientificname,
# value.var = "occurrence", fill = 0
# )

# # Replace NA with 0 for non-detections
# data_wide[is.na(data_wide)] <- 0

# Let's try something different
tmp <- matrix(0, nrow = nrow(data_all2), ncol = length(all_brbs))
rownames(tmp) <- data_all2$sampling_event_identifier
colnames(tmp) <- all_brbs

pb <- txtProgressBar(min = 0, max = 100, style = 3)
for (i in 1:nrow(data_all2)) {
    event <- data_all2$sampling_event_identifier[i]
    species <- data_all2$scientificname[i]
    tmp[event, species] <- 1
    setTxtProgressBar(pb, (i / nrow(data_all2)) * 100)
}
t2 <- Sys.time()
# This might be expensive - don't need this for the moment
tmp2 <- aggregate(tmp, by = list(rownames(tmp)), FUN = sum)
t3 <- Sys.time()
# Save the data
ebird_final <- tmp2
print(num)
save(ebird_final, file = paste0("~/Downloads/ebird_transform", num, ".Rdata"))

# Compile

n <- 10
store <- list()
for (i in 1:10) {
    print(i)
    load(paste0("~/Downloads/ebird_transform", i, ".Rdata"))
    store[[i]] <- ebird_final
}

data <- do.call(rbind, store)

# Now we need to match the effort/lat-lon information to the ebird data
tmp <- data_all[match(data$Group.1, data_all$sampling_event_identifier), c("sampling_event_identifier", "latitude", "longitude", "year", "number_observers", "effort_distance_km", "duration_minutes")]

data_final <- cbind(tmp, data)

# Remove the group.1 column
data_final$Group.1 <- NULL

raw <- data_final

# Save
save(raw, file = "~/phylo-sdms/phyloproj/raw_data/hummingbirdSA/ebird_raw-2024.Rdata")
