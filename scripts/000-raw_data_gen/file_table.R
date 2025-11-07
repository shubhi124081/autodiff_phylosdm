root <- "/Users/ss4224/phylo-sdms2"
raw_directory <- file.path(root, "raw_data/hummingbirdSA")
CLUSTER <- "Coeligena"
files <- dir(file.path(raw_directory, CLUSTER), pattern = "*\\.Rdata$", full.names = TRUE)
files <- substr(files, nchar(raw_directory) + nchar(CLUSTER) + 3, nchar(files))
files <- sub("\\.Rdata$", "", files)
files <- files[grepl("^y\\d+", files)]

# Split
bits <- strsplit(files, "_")
df <- data.frame(
    def_lev = sapply(bits, function(x) x[1]),
    species = sapply(bits, function(x) paste(x[2], x[3], sep = "_")),
    stringsAsFactors = FALSE
)
df <- unique(df)
df <- df[order(df$species), ]

write.csv(df, "~/Downloads/possible_art_exps.csv", row.names = FALSE)
