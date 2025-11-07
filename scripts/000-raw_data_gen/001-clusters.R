# Libraries & functions ----
source("~/phylo-sdms/phyloproj/scripts/000-raw_data_gen/010-raw_data_functions.R") # nolint
source("~/phylo-sdms/phyloproj/scripts/000-phyloGenie_functions.R")

# Function to match species names
match_species_names <- function(raw, tree, start_col) {
  ebird_species <- colnames(raw[, start_col:ncol(raw)])
  tree_species <- tree$tip.label

  matched_species <- intersect(ebird_species, tree_species)
  unmatched_species <- setdiff(ebird_species, tree_species)

  list(matched = matched_species, unmatched = unmatched_species)
}

# Make clusters function
makeSpCluster <- function(GENERA, DIST, CUTOFF) {
  genera2 <- GENERA
  spList <- list()
  count <- 1
  nn <- ""

  for (i in seq_along(GENERA)) {
    genus <- GENERA[i]
    if (genus %in% genera2) {
      ind <- grep(genus, rownames(DIST))

      cutoff <- if (genus %in% c("Androdon", "Eutoxeres")) 1.5 else CUTOFF
      sps <- spByDist(DIST, cutoff, ind[1])
      spList[[count]] <- sps
      nn <- append(nn, genus)

      inc <- getGenus(sps)
      print(genus)
      leftBehindWarning(rownames(DIST), sps)

      rmv <- which(genera2 %in% inc)
      genera2 <- genera2[-rmv]
      count <- count + 1
    }
  }
  nn <- nn[-1]
  names(spList) <- nn
  return(spList)
}

# Set-up ----
dpath <- "~/phylo-sdms2/raw_data/hummingbirdSA/"
dataset <- "ebird_raw-2024.Rdata"
phy <- "mcguire_tree.Rdata"

# Load data sets ----
contents <- load(file = paste0(dpath, dataset)) # loads obj raw
colnames(raw) <- gsub(" ", "_", colnames(raw))
print(paste0("Raw eBird dataset has: ", ncol(raw)))

# Load tree
contents <- load(file = paste0(dpath, phy)) # loads obj tree
print(paste0("Raw phylogeny has: ", length(tree$tip.label)))

# Another dose of harmonization
# Load the harmonization data
recon <- read.csv(file = paste0(dpath, "ebird_mctree_harmonization.csv"))
tree$tip.label[order(tree$tip.label)] <- recon$mc_tree_edited

# Match species names -----
# Determine the starting column for species names
species_start_col <- which(sapply(colnames(raw), function(x) grepl("^[A-Z][a-z]+_[a-z]+$", x)))[1]

# Harmonize species names between raw data and tree
species_match <- match_species_names(raw, tree, species_start_col)
print(paste0("Matched species: ", length(species_match$matched)))
print(paste0("Unmatched species: ", length(species_match$unmatched)))

# Normalize branch lengths
tree$edge.length <- tree$edge.length / max(phytools::nodeHeights(tree)[, 2]) * 1

# Harmonize
y <- raw[, species_start_col:ncol(raw)]
y <- y[, species_match$matched]
tree <- ape::keep.tip(tree, species_match$matched)

# Complete cases only
cc <- which(complete.cases(raw) == TRUE)
raw <- raw[cc, ]

# Create cood df
cood <- raw[, c("latitude", "longitude")]

# Configuration ----
MAX_SP <- 20
CUTOFF <- 1

# Get genera and distance ----
genera <- getGenus(species_match$matched)
genera <- genera[order(genera)]
dist <- ape::cophenetic.phylo(tree)

# Make clusters round 1 -----
spList <- makeSpCluster(genera, dist, CUTOFF)
nn <- names(spList)

# Rework clusters -----
# Pick clusters that are greater than 20 species

for (i in seq_along(spList)) {
  len <- length(spList[[i]])

  # Check if length of any spList elements exceeds 20 (reasonable computational limit)
  if (len > 20) {
    sps <- spList[[i]]
    sps_name <- names(spList)[i]

    # Determine split factor (2 or 3) based on length
    split_factor <- ifelse(len / 2 > MAX_SP, 3, 2)

    # Generate new names for split elements
    new_names <- paste0(sps_name, 1:split_factor)

    # Compute split indices
    cuts <- floor(seq(1, len, length.out = split_factor + 1))

    # Update original list and create new elements
    spList[[i]] <- sps[cuts[1]:cuts[2]]
    for (j in 2:split_factor) {
      spList[[new_names[j]]] <- sps[(cuts[j] + 1):cuts[j + 1]]
    }
  }
}

# Save -----
save(spList, file = file.path(dpath, "spList.Rdata"))
save(tree, file = paste0(dpath, phy))

# Update - Archilochus2 does not run with the current setup
# I'm manually splitting it in Archilochus2, 3, and 4
# But then I reversed it...

contents <- load(file = file.path(dpath, "spList.Rdata")) # loads obj spList
spList$Archilochus2 <- c(spList$Archilochus2, spList$Selasphorus)
spList$Selasphorus <- NULL
save(spList, file = file.path(dpath, "spList.Rdata"))

# Archive

# groups <- list(
#   "Adelomyia" = c(),
#   "Archilochus" = c(4, 18, 27),
#   "Phaethornis" = c(15, 25)
# )
# w <- which(names(spList) %in% names(groups))
# # SpList corrections
# correctSpCluster <- function(SPLIST, GROUPS, INDEX) {
#   # Set up
#   SPLIST <- spList
#   GROUPS <- groups
#   INDEX <- w
#   count <- length(SPLIST) + 1

#   # First loop over index where the spList needs to be broken down
#   for (index in INDEX) {
#     sps <- SPLIST[[index]]
#     splits <- GROUPS[[names(SPLIST)[index]]]
#     # Second loop over the length of the splits with each cluster
#     for (j in seq_len(length(splits))) {
#       if (j == 1) {
#         SPLIST[[count]] <- sps[1:splits[j]]
#       } else {
#         SPLIST[[count]] <- sps[(splits[j - 1] + 1):splits[j]]
#       }
#       count <- count + 1
#     }
#   }
#   return(SPLIST)
# }

# spList <- correctSpCluster(spList, groups, w)
# spList <- spList[-w]
# nn <- nn[-w]
# nn <- c(
#   nn, "Opisthoprora", "Metallura", "Adelomyia", "Amazilia",
#   "Eupherusa", "Archilochus", "Chaetocercus", "Selasphorus", "Phaethornis",
#   "Phaethornis2"
# )
# names(spList) <- nn

# # Remove NAs from spList
# spList <- lapply(spList, function(x) x[!is.na(x)])
# # Nope, we lost like 20 species?

# for (i in seq_len(length(spList))) {
#   len <- length(spList[[i]])
#   # Check if length of any spList elements exceed 20 (this is reasonable computational limit)
#   if (len > 20) {
#     sps <- spList[[i]]
#     sps_name <- names(spList)[i]
#     split <- length(sps) / 2
#     # If split by 2 is not enough, split by 3
#     if (split > 20) {
#       # New names
#       sps_newname1 <- paste0(sps_name, "2")
#       sps_newname2 <- paste0(sps_name, "3")
#       # Length cut off is length of spList[element] divided by 3
#       cutoff <- length(sps) / 3
#       # Round down
#       cuts <- floor(seq(1, length(sps), length.out = 4))
#       # Cut original
#       spList[[i]] <- sps[cuts[1]:cuts[2]]
#       # Make new dummy list elements with rest of the cuts
#       spList$XX <- sps[(cuts[2] + 1):cuts[3]]
#       spList$XX2 <- sps[(cuts[3] + 1):cuts[4]]
#       # Change name of dummy list elements
#       chgn <- grep("XX", names(spList))
#       names(spList)[chgn] <- c(sps_newname1, sps_newname2)
#     } else {
#       # New names
#       sps_newname <- paste0(sps_name, "2")
#       # Length cut off is length of spList[element] divided by 2
#       cutouff <- length(sps) / 2
#       # Round down
#       cuts <- floor(seq(1, length(sps), length.out = 3))
#       # Cut original
#       spList[[i]] <- sps[cuts[1]:cuts[2]]
#       # Make dummy list elements
#       spList$XX <- sps[(cuts[2] + 1):cuts[3]]
#       # Change name of dummy list elements
#       chgn <- grep("XX", names(spList))
#       names(spList)[chgn] <- c(sps_newname)
#     }
#   }
# }
