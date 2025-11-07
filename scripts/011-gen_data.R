##################### Set-up #####################
# Load all libraries at once
rm(list = ls())

# Some libs
library(ape)
library(phytools)

# UGH
root <- "~/phylo-sdms2"

# Filepaths & scripts
script_directory <- file.path(root, "scripts")
data_directory <- file.path(root, "data")
res_directory <- file.path(root, "res")
job_directory <- file.path(root, "jobs")
raw_directory <- file.path(root, "raw_data")

# All the functions we need
source(file.path(script_directory, "000-phyloGenie_functions.R"))

##################### Read Files In  #####################
# Source simulation_build_functions.R
# source(file.path(script_directory, "011-sim_build_fns.R"))

# Configuration parameters for data generation
config <- yaml::read_yaml(file.path(script_directory, "000-config.yaml"),
  as.named.list = TRUE
)
data_list <- config$data
data_sim <- config$data_sim
# shortname is going to be the name of the experiment folder and the shorthand for
# file IDs
shortname <- paste0(data_list$exp_root, "_", data_list$exp_id, "_", data_list$cluster)

##################### Commit config file now #####################

if (data_list$autocommit) {
  # RStudio wd has to be set to scripts dir first
  setwd(script_directory)
  system(paste0("git add ", script_directory, "/000-config.yaml"))
  system(paste0(
    "git commit -m '(auto-commit) config params for ",
    shortname, "'"
  ))
  system("git push")
}

##################### If sim == TRUE  #####################

if (data_sim$sim) {
  config$phylo$mean_root_vector <- yamlConvert(vec = TRUE)
  config$phylo$covar_traits <- yamlConvert(vec = FALSE)

  phylo <- config$phylo
  mod <- config$mod
  nseeds <- data_sim$nseeds

  ##################### Loops #####################
  parent <- data_sim$var$parent
  features <- data_sim$var$feature
  # Safety check : make sure toggle features match mod inputs
  if (!all(features %in% names(config[[parent]]))) {
    warning("Your selected features don't match any mod obj inputs")
  }

  # Make parameter table
  d1 <- makeParameterTbl()

  #### Read in a phylo ####
  # Loads object named tree
  # NOTE - when run interactive, filepath returned from getwd()
  # root <- gsub("scripts", "", getwd())
  load(file.path(root, data_sim$tree_relpath))
  tree$tip.label <- LETTERS[seq_len(length(tree$tip.label))]

  for (i in seq_len(nrow(d1))) {
    print(paste0("VAR: ", colnames(d1), " VALUE: ", d1[i, ]))
    # Preserve mod by saving mod2 obj
    parent2 <- config[[parent]]
    ind <- which(names(parent2) %in% features)

    # Change values of selected features
    for (ii in seq_len(length(ind))) {
      parent2[[ind[ii]]]$value <- d1[i, names(parent2)[ind[ii]]]
    }

    # Reps
    for (seed in seq_len(nseeds)) {
      print(paste0("SEED NO: ", seed))

      # Make dataframe based on changed values
      if (parent == "phylo") {
        params <- parent2
      } else {
        params <- phylo
      }

      data <- simulateAbundanceData1(TREE = tree, simulation_params = params)
      grid <- createGrid(
        user_grid = data_sim$user_grid,
        nsamples = data_sim$num_rec,
        path = file.path(root, data_sim$grid_relpath)
      )

      # Simulates an OU/BM process
      all_data <- simulateAbundanceData_new(
        grid = grid,
        nsamples = data_sim$num_rec,
        binary = data_sim$binary,
        print = TRUE,
        tree,
        vars = c("CHELSA_bio_1.tif", "CHELSA_bio_13.tif"),
        ntraits = params$ntraits,
        sigma2 = params$sigma2$value,
        alpha = params$alpha$value,
        epath = "~/env"
      )

      # Organize data into long form instead of wide form
      Ytrue <- all_data$Ytrue
      X <- all_data$X
      true_betas <- all_data$true_betas
      train_index <- all_data$train_index

      # If binary, need to convert everything in this directory to presence-only
      # But keep binary for eval
      if (data_sim$binary) {
        sp <- colnames(Ytrue)
        sp <- sp[-which(sp == "rwid")]
        Ytrue_vector <- unlist(lapply(sp, function(species) {
          rep(species, sum(Ytrue[, species] == 1))
        }))
        row_ids <- unlist(lapply(sp, function(species) {
          which(Ytrue[, species] == 1)
        }))
        Ytrue_pp <- data.frame(species = Ytrue_vector, row_id = row_ids)

        X <- X[row_ids, ]
        Ytrain <- Ytrue_pp[which(Ytrue_pp$row_id %in% train_index), ]
        Ytest <- Ytrue_pp[-which(Ytrue_pp$row_id %in% train_index), ]
        X <- X[, -c(which(colnames(X) == "rwid"))]
        Xtrain <- X[train_index, ]
        Xtest <- X[-train_index, ]
        Ytrue <- Ytrue_pp
      } else {
        # If continuous
        # Subsetting matrices into test and train
        Ytrain <- all_data$Ytrain
        Ytest <- all_data$Ytest
        Ytrain <- Ytrain[, -c(which(colnames(Ytrain) %in% c("rwid", "tot")))]
        Ytest <- Ytest[, -c(which(colnames(Ytest) %in% c("rwid", "tot")))]
        Y_occ <- Ytrain
        Y_occ[Y_occ < 0] <- 0
        Xtrain <- all_data$Xtrain
        Xtrain <- Xtrain[, -c(which(colnames(Xtrain) == "rwid"))]
        Xtest <- all_data$Xtest
        Xtest <- Xtest[, -c(which(colnames(Xtest) == "rwid"))]
        true_betas <- all_data$true_betas
        train_index <- all_data$train_index
        X <- X[, -c(which(colnames(X) == "rwid"))]
        Ytrue <- all_data$Ytrue[, -c(which(colnames(all_data$Ytrue) == "rwid"))]
      }
      true_tree <- tree

      # Making matrices for stan code
      ns <- length(unique(Ytrue$species))
      nsamples_train <- nrow(Ytrain)
      L <- matrix(ape::vcv(tree, corr = TRUE), ns, ns)
      distance_matrix <- ape::cophenetic.phylo(tree)
      I <- matrix(0, ns, ns)
      diag(I) <- 1
      I_y <- matrix(0, nsamples_train, nsamples_train)
      diag(I_y) <- 1

      # Putting all data together
      print(paste0("NO. OF SPECIES = ", ns))
      standata <- list(
        N = nsamples_train,
        K = ncol(Xtrain),
        J = ns,
        E = 1,
        N_occ = nsamples_train,
        y = as.matrix(Ytrain),
        y_occ = as.matrix(Ytrain),
        y_test = as.matrix(Ytest),
        x_test = as.matrix(Xtest),
        x = Xtrain,
        x_occ = Xtrain,
        x_tilde = Xtrain,
        N_tilde = nrow(Xtrain),
        L_dist = distance_matrix,
        L_corr = L,
        beta_mu = rep(0, ncol(Ytrain)),
        true_betas = true_betas,
        I = I,
        Iy = I_y,
        alpha = 1,
        st_devs = rep(1, ncol(Ytrain)),
        U = mod$U$value,
        a = mod$a$value,
        b = mod$b$value,
        tr = tree,
        tiptraits = data$tip_traits,
        Y_true = Ytrue,
        X_true = X,
        train_index = train_index
      )

      if (parent == "mod") {
        standata[!is.na(match(names(standata), features))] <-
          sapply(parent2, function(x) {
            return(x$value)
          })
      }

      # Setting up file name
      mid <- paste0(data_sim$var$feature, collapse = "_")
      mid_vals <- paste0(round(d1[i, ], 2), collapse = "_")

      # Create a data directory and filepath
      data_dir <- file.path(data_directory, shortname)

      if (!dir.exists(data_dir)) {
        dir.create(data_dir)
      }

      filename <- paste0(
        data_sim$experiment_name, "-var_", mid, "_",
        mid_vals, "_rep_", seed, ".Rdata"
      )

      # Adding config params to the experiment file
      config_save <- config
      config_save[[parent]] <- parent2
      everything <- list(
        data = standata,
        config = config_save
      )
      # Save
      save(everything, file = paste0(
        data_dir, "/",
        filename
      ))
    }
  }
}

##################### If sim == FALSE  #####################

if (!data_sim$sim) {
  mod <- config$mod
  fsp <- data_list$focal_sp
  def_lev <- data_list$def_level
  dataset_name <- data_list$raw_data
  cluster <- data_list$cluster
  wd <- file.path(raw_directory, dataset_name, cluster)
  repno <- data_list$repno
  nrep <- data_list$nrep
  offset_vars <- data_list$offsets
  if (nrep > 1) {
    rep_sequence <- seq_len(nrep)
  } else {
    rep_sequence <- repno
  }

  # Y <- read.csv(paste0(wd, "/y.csv"))
  # X <- read.csv(paste0(wd, "/x.csv"))
  # loads an obj called store
  # TODO: Add fold or rep here!

  for (jj in seq_len(length(def_lev))) {
    def_lev1 <- def_lev[jj]

    # Load the full dataset
    contents <- load(file.path(wd, paste0(def_lev1, "_", fsp, "_run_files.Rdata")))

    # Make sure subsetting and NA removal doesn't mess everything up
    # Three matricies
    # x and cood are site level, y is species-site level
    store_y <- store$y
    store_x <- store$x
    store_cood <- store$cood

    for (i in rep_sequence) {
      # Load the train-test indices
      contents <- load(file.path(wd, paste0(def_lev1, "_", fsp, "_indices.Rdata")))

      idx <- list_of_indices[[i]]
      idx_train <- idx$training_sites
      idx_train <- idx_train[order(idx_train)]

      # Subset y
      y <- store_y[store_y$site %in% idx_train, ]
      # Match site IDs to rows in x and cood
      x <- store_x[y$site, , drop = FALSE]
      cood <- store_cood[y$site, , drop = FALSE]

      # loads an obj called tree
      # contents <- load(paste0(wd, "/tree.Rdata"))
      tree <- store$tree

      if (class(tree) == "phylo") {
        distance_matrix <- ape::cophenetic.phylo(tree)
        # L <- matrix(ape::vcv(tree, corr = T), ncol(Y),ncol(Y))
      } else {
        distance_matrix <- tree
      }

      # Separate y matricies
      y_species_idx <- y[, "species"]
      y_values <- y[, "count"]
      y_sites <- y[, "site"]
      if (data_list$thin) {
        y_values <- ifelse(y_values > 0, 1, 0)
      }
      if (min(y_values) != 0) {
        stop("y_values should be zero or one for presence-absence data.")
      }

      # Construct offset matrix
      # offset_vars <- c("duration", "distance", "num_observers")
      # if ("soft_clip" %in% colnames(x)) {
      #   offset_vars <- c(offset_vars, "soft_clip")
      # }

      # Ensure all offset_vars exist in x
      missing_vars <- setdiff(offset_vars, colnames(x))
      if (length(missing_vars) > 0) {
        stop(paste("Missing offset variables in x:", paste(missing_vars, collapse = ", ")))
      }
      offset <- x[, offset_vars, drop = FALSE]
      offset <- cbind(offset, site = y$site)
      # Remove only those offset columns that actually exist in x
      cols_to_remove <- intersect(colnames(x), offset_vars)
      if (length(cols_to_remove) > 0) {
        x <- x[, !(colnames(x) %in% cols_to_remove), drop = FALSE]
      }
      colnames(offset) <- c(offset_vars, "site")

      # Calculate offset as product of duration, distance, num_observers (and soft_clip if present)
      offset_matrix <- offset
      offset_vars_to_mult <- setdiff(colnames(offset_matrix), "site")
      if (length(offset_vars_to_mult) == 1) {
        offset_product <- offset_matrix[, offset_vars_to_mult]
      } else {
        offset_product <- apply(offset_matrix[, offset_vars_to_mult, drop = FALSE], 1, prod)
      }
      offset_multiplied <- data.frame(offset = offset_product, site = offset_matrix[, "site"])
      offset_vector <- offset_multiplied[match(y_sites, offset_multiplied$site), "offset"]
      offset_vector[is.na(offset_vector)] <- 1
      # Offset vector is logged later in 023-main.R

      # Get ready for stan
      x_matrix <- x

      # Create stan data
      stan_data <- list(
        N = nrow(x_matrix),
        J = length(unique(y_species_idx)),
        K = ncol(x_matrix),
        N_obs = length(y_values),
        species = y_species_idx,
        X = x_matrix,
        y = as.integer(y_values),
        D_phylo = distance_matrix,
        offset = offset_vector,
        tr = tree
      )

      # Create a data directory and filepath
      data_dir <- file.path(data_directory, shortname)

      if (!dir.exists(data_dir)) {
        dir.create(data_dir)
      }

      everything <- list(
        data = stan_data,
        config = config
      )
      # Naming structure -
      # EXPROOT_EXPID_CLUSTER_FSP_DEFLEV_REPNO.Rdata
      # Experiment root is "happy", "grumpy" etc.
      # Experiment ID is some sort of note about this round e.g. LGCP1 (first round of LGCP)
      # Cluster is the cluster of species
      # Def level is the level of deficiency
      # Model is the model used
      # Repno is the rep number

      filename <- paste0(shortname, "_", fsp, "_", def_lev1, "_", i, ".Rdata")
      save(everything, file = paste0(
        data_dir, "/",
        filename
      ))
    }
  }
}
system(paste0("echo ", paste0(shortname)))
