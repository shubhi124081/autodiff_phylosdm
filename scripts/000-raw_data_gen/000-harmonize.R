# Script to harmonize the raw data from the different sources
# Specifically for the hummingbird data from eBird and GBIF
# and the species list from the McGuire phylogeny

# Mcguire phylogeny
mc_tree <- ape::read.tree("~/Downloads/hum294_correctedIDs.tre")
# loads obj tree
contents <- load("~/phylo-sdms/phyloproj/raw_data/hummingbirdSAtree.Rdata")

# Need to simplify species names from mc_tree
mc_sps <- mc_tree$tip.label
# Remove everything after the second "."
mc_sps <- sub("^(\\w+\\.\\w+)\\..*$", "\\1", mc_sps)
mc_sps <- sub("\\.", "_", mc_sps) # Replace the first "." with "_"

# Oh no - it looks like there is signficiant cleaning to be done with the
# species names from the McGuire phylogeny

# First, let's make sure all the species names have the first character
# capitalized and the rest not
mc_sps <- tolower(mc_sps)
mc_sps <- sub("^(\\w)", "\\U\\1", mc_sps, perl = TRUE)

# Now, some of the genus names are shortened
# We need to match them to the closest genus name in tree$tip.label and replace
# Extract genus names from tree$tip.label
tree_genera <- unique(sub("_.*", "", tree$tip.label))

# Function to find the closest match by checking if the genus is a prefix
closest_genus <- function(genus, genera_list) {
    matches <- genera_list[startsWith(genera_list, genus)]
    if (length(matches) > 0) {
        return(matches[1])
    } else {
        return(NA)
    }
}

# Match genus names in mc_sps to the closest genus name from tree$tip.label
mc_genera <- sub("_.*", "", mc_sps)
closest_genera <- sapply(mc_genera, closest_genus, genera_list = tree_genera)

# Read in the spreadsheet that does the harmonization
recon <- read.csv("~/phylo-sdms/phyloproj/raw_data/hummingbirdSA/taxanomic_harmonization.csv")
recon$New_name <- sub(" ", "_", recon$New_name)

# Match mc_sps_raw and New_name
# Change the names in mc_sps to the New_name counterpart
# Then will be ready to create new spList and merge with ebird data
mm <- match(mc_sps, recon$Mc_sps_raw)
mc_sps[!is.na(mm)] <- recon$New_name[mm[!is.na(mm)]]

# Change the name of the tree tip labels
mc_tree$tip.label <- mc_sps
tree <- mc_tree
save(tree, file = "~/phylo-sdms/phyloproj/raw_data/hummingbirdSA/mcguire_tree.Rdata")

# Archive

# # Create a dataframe with mc_sps and the closest genus name from tree$tip.label
# harmonized_df <- data.frame(mc_sps = mc_sps, closest_genus = closest_genera)

# # Replace mc_sps genus with closest_genus where closest_genus is not NA
# harmonized_df$mc_sps <- ifelse(!is.na(harmonized_df$closest_genus),
#     paste(harmonized_df$closest_genus,
#         sub(".*_", "", harmonized_df$mc_sps),
#         sep = "_"
#     ),
#     harmonized_df$mc_sps
# )


# # Let's look at which tree$tip.label species are not in the harmonized_df
# missing_sps <- tree$tip.label[!tree$tip.label %in% harmonized_df$mc_sps]

# # Can we look at the closest strings to the missing species?
# # Function to find the closest string match
# closest_match <- function(string, string_list) {
#     distances <- adist(string, string_list)
#     closest_index <- which.min(distances)
#     return(string_list[closest_index])
# }

# # Create a dataframe with missing_sps and their closest string match from harmonized_df$mc_sps
# closest_matches <- sapply(missing_sps, closest_match, string_list = harmonized_df$mc_sps)
# missing_df <- data.frame(missing_sps = missing_sps, closest_match = closest_matches)

# # View the dataframe
# print(missing_df)

# # Okay, let's read in the master hummingbird species list
# master_sps <- read.csv("~/phylo-sdms/phyloproj/raw_data/IOC_hummingbird_master_list.csv")
# master_sps$sp_name <- paste0(master_sps$Genus, "_", master_sps$Species)
# # Alphabtical order
# master_sps <- master_sps[order(master_sps$sp_name), ]
# # Alphabetize mc_sps
# mc_sps <- sort(mc_sps)
