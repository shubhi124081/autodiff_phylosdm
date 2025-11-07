path <- "~/phylo-sdms2/raw_data/hummingbirdSA"
contents <- load("~/phylo-sdms2/raw_data/hummingbirdSA/spList.Rdata") # loads obj spList
vec_of_names <- c("t2003_t2006", "t2003_t2009", "t2003_t2012", "t2003_e2015", "t2003_t2018", "t2019_t2024")

for (sp in names(spList)) {
    sp_dir <- file.path(path, sp)
    if (dir.exists(sp_dir)) {
        for (name in vec_of_names) {
            old_file <- file.path(sp_dir, paste0(name, "_ALL_run_files.Rdata"))
            if (file.exists(old_file)) {
                # Change e2003-e2006 -> t2003_t2006
                new_name <- sub("^(t\\d{4})_(t\\d{4})", "\\1-\\2", name)
                new_file <- file.path(sp_dir, paste0(new_name, "_ALL_run_files.Rdata"))
                file.rename(old_file, new_file)
            }
        }
    }
}
