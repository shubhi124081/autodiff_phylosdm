# Raw data wrangling functions
# Organize data as cood-y

organizeY <- function(dpath, name) {
  if (name == "hummingbird_sa") {
    print("1. Loading data...")
    raw <- read.csv(file.path(dpath, paste0(name, ".csv")))

    print("2. Organizing data...")
    # Replace "." in sp names with "_"
    # Species start at column 8
    nm <- colnames(raw)
    st <- 8
    colnames(raw)[st:ncol(raw)] <- gsub("\\.", "_", nm[st:length(nm)])
    spnm <- colnames(raw)[st:ncol(raw)]

    # Grab coordinates and change "x", "y" to "lon", "lat" respectively
    cood <- which(nm %in% c("y", "x"))
    colnames(raw)[cood] <- c("lat", "lon")

    return(raw[, c("lat", "lon", "observation_date", "num.observers", "duration", "distance", spnm)]) # nolint
  }
}

# Check phylo sp and occ sp are the same

equalizeSP <- function(tree, y) {
  # y needs to be just species records
  tsp <- tree$tip.label
  ysp <- colnames(y)

  # Which y sp NOT IN t sp
  ky <- which(!(ysp %in% tsp))
  if (!identical(ky, integer(0))) y <- y[, -c(ky)]

  # Which t sp NOT IN y sp
  ysp <- colnames(y)
  kt <- which(!(tsp %in% ysp))
  if (!identical(kt, integer(0))) tree <- ape::drop.tip(tree, kt)

  # Final check
  tsp <- tree$tip.label
  if (all(ysp %in% tsp) & all(tsp %in% ysp)) {
    print("All species names match")
  } else {
    print("Something went wrong - debug with y <- Y; tree <- TREE")
  }

  return(list(y = y, tree = tree))
}


# X covariates

organizeCovar <- function(epath, env.files, env.crs, cood) {
  env.paths <- file.path(epath, env.files)
  colnames(cood) <- c("y", "x")
  cood <- cood[, c("x", "y")]
  extent <- raster::extent(cood)

  env <- raster::stack()
  for (env.path in env.paths) {
    print(paste0("Loading ", env.path, " ..."))
    env.layer <- raster::raster(env.path)

    # If needed, reproject the env layer
    if (!raster::compareCRS(raster::crs(env.layer), env.crs)) {
      print(paste0(" Reprojecting..."))
      env.layer <- raster::projectRaster(env.layer, env.crs)
    }

    if (!is.null(extent)) {
      env.crop <- raster::crop(env.layer, extent)
    }
    env <- raster::addLayer(env, env.crop)
  }

  cood.mat <- do.call(cbind, cood)

  annotations <- raster::extract(env, cood.mat,
    method = "bilinear"
  )
  return(annotations)
}

# Function for creating a reference coord set

makeRefCoods <- function(REF_RAST) {
  REF_RAST[] <- 1
  area <- rnaturalearth::ne_countries(
    scale = "medium",
    returnclass = "sf", continent = c("north america", "south america")
  )
  area_vect <- terra::vect(area)
  # Mask
  REF_RAST <- terra::mask(REF_RAST, area_vect)
  w <- terra::as.data.frame(REF_RAST, xy = TRUE)
  # Return reference set of coordinates
  return(w)
}

gimmeRatios <- function(YFSP, YCOOD, REF_RAST) {
  # Parition the YFSP
  # First quick check
  if (length(YFSP) != nrow(YCOOD)) {
    error("Something is wrong - Y and COOD are not equal")
  }
  ind1 <- which(YFSP == 1)
  ind0 <- which(YFSP == 0)

  # Dealing with absences
  # Create empty raster for recording abscenes
  l <- dim(REF_RAST)[1] * dim(REF_RAST)[2]
  er0 <- REF_RAST
  er0[] <- 0
  # Count how many absence coordinates fall into each cell
  counts0 <- table(terra::cellFromXY(er0, YCOOD[ind0, ]))
  c0 <- data.frame("n" = as.numeric(names(counts0)), "count" = counts0)
  df0 <- data.frame("n" = 1:l, "count" = rep(0, l))
  df0[match(c0$n, df$n), "count"] <- c0$count.Freq
  # Annotate the empty raster with the counts
  er0[] <- df0$count
  tmp0 <- terra::as.data.frame(er0, xy = TRUE)

  # Dealing with presences
  # Create empty raster for recording presences
  er1 <- REF_RAST
  er1[] <- 0
  counts1 <- table(terra::cellFromXY(er1, YCOOD[ind1, ]))
  c1 <- data.frame("n" = as.numeric(names(counts1)), "count" = counts1)
  df1 <- data.frame("n" = 1:l, "count" = rep(0, l))
  df1[match(c1$n, df1$n), "count"] <- c1$count.Freq
  er1[] <- df1$count
  tmp1 <- terra::as.data.frame(er1, xy = TRUE)

  err <- er1 / er0
  # tmpR <- terra::as.data.frame(err, xy = TRUE)
  return(err)
}

# Functions
createSampleRast <- function(EX, AREA_VECT, EPATH, ENV_FILE) {
  env_path <- file.path(EPATH, ENV_FILE)
  env <- terra::rast(env_path)
  env_crop <- terra::crop(env, ex)
  r <- terra::mask(env_crop, AREA_VECT)

  er <- terra::rast(terra::ext(r), resolution = terra::res(r))
  terra::crs(er) <- terra::crs(r) # sample raster

  return(er)
}


createEffortRast <- function(SAMPLE_RAST, EC, VAR,
                             TIME = NULL, PATH = NULL) {
  if (!is.null(TIME)) {
    EC <- EC[which(EC$year <= TIME), ]
  }
  e1_vect <- terra::vect(EC[, c("lat", "lon", VAR)],
    geom = c("lon", "lat")
  )
  e1 <- terra::rasterize(e1_vect, SAMPLE_RAST, field = VAR, fun = sum)

  if (!is.null(PATH)) {
    if (!is.null(TIME)) {
      name <- paste0("effort_rast_", VAR, "_e", TIME, ".tiff")
    } else {
      name <- paste0("effort_rast_", VAR, "_FULL.tiff")
    }
    terra::writeRaster(e1, file = file.path(PATH, name))
  }
  return(e1)
}


# Functions

annotateCoods <- function(EPATH, ENV_FILES, ENV_CRS, COOD) {
  env_paths <- file.path(EPATH, ENV_FILES)
  cood <- COOD[, c("lon", "lat")]
  colnames(cood) <- c("x", "y")
  extent <- terra::ext(getExtent(cood))

  env <- list()
  for (i in seq_len(length(ENV_FILES))) {
    env_path <- env_paths[i]
    print(paste0("Loading ", env_path, " ..."))
    env_layer <- terra::rast(env_path)

    # Reproject
    print(paste0(" Reprojecting..."))
    env_layer <- terra::project(env_layer, ENV_CRS)

    if (!is.null(extent)) {
      env_crop <- terra::crop(env_layer, extent)
    }
    env[[i]] <- env_crop
  }
  cood_mat <- do.call(cbind, cood)

  ###### 6. Annotate coordinates with envs #######
  a <- list()
  for (i in seq_len(length(ENV_FILES))) {
    print(i)
    a[[i]] <- terra::extract(x = env[[i]], y = cood_mat)
  }

  annotations <- do.call(cbind, a)
  if (any(grepl("\\.tif", ENV_FILES))) {
    if (any(grepl("\\.tiff", ENV_FILES))) {
      xx <- gsub("\\.tiff", "", ENV_FILES)
    } else {
      xx <- gsub("\\.tif", "", ENV_FILES)
    }
  }
  colnames(annotations) <- xx
  return(annotations)
}

# Function to take in a data frame of coordinates and crop by
# an extent, returning a cropped data frame of coordinates

makeOccurrenceMap <- function(YC, FSP) {
  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  if (all(c("lat", "lon") %in% colnames(YC))) {
    ind_lat <- which(colnames(YC) == "lat")
    ind_lon <- which(colnames(YC) == "lon")
    colnames(YC)[c(ind_lon, ind_lat)] <- c("x", "y")
  }
  COOD <- YC[, c("x", "y")]
  ex <- raster::extent(COOD)
  xmin <- ex[1]
  xmax <- ex[2]
  ymin <- ex[3]
  ymax <- ex[4]
  df <- YC[, c("x", "y", FSP)]
  df$col <- ifelse(df[, 3] == 1, "Presence", "Absence")
  ggplot() +
    geom_point(data = df[which(df[, 3] == 0), ], aes(
      x = x,
      y = y
    ), shape = 21, size = 3, fill = "#D98949") +
    geom_point(
      data = df[which(df[
        ,
        3
      ] == 1), ], aes(x = x, y = y), shape = 21, size = 3,
      fill = "#7CB781"
    ) +
    geom_sf(
      data = world, fill = NA,
      col = "black"
    ) +
    coord_sf(xlim = c(xmin, xmax), ylim = c(
      ymin,
      ymax
    )) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(fill = "") +
    ggtitle(paste0(
      "Occurrence map for ",
      FSP
    ), subtitle = paste0(
      "Number of presence points = ",
      sum(df[, 3]), "\nNumber of Absence points = ",
      (nrow(df) - sum(df[, 3]))
    )) +
    ylab("Latitude") +
    xlab("Longitude")
}

cropDFCood <- function(YCOOD, EX, NAMES = c("lon", "lat")) {
  if (class(EX) != "SpatExtent") EX <- terra::ext(EX)

  cood_vect <- terra::vect(YCOOD, geom = c("lon", "lat"))
  cood_crop <- terra::crop(cood_vect, EX)
  y_crop <- terra::values(cood_crop)
  cood_crop <- terra::geom(cood_crop)[, c("x", "y")]
  yc <- cbind(cood_crop, y_crop)
  return(yc)
}

duplicatedCood <- function(DF) {
  DF$str <- paste0(DF$x, "_", DF$y)
  ind <- which(duplicated(DF$str))
  return(ind)
}

# This is gimmeRatios function
gimmeRatios <- function(YFSP, YCOOD, REF_RAST, NAMES = c("x", "y")) {
  if (length(YFSP) != nrow(YCOOD)) {
    error("Something is wrong - Y and COOD are not equal dimensions")
  }
  ind1 <- which(YFSP == 1)
  ind0 <- which(YFSP == 0)
  c0 <- YCOOD[ind0, ]
  c0_vect <- terra::vect(c0[, NAMES], geom = NAMES)
  er0 <- terra::rasterize(c0_vect, REF_RAST, fun = length)
  er0[is.na(er0)] <- 0
  c1 <- YCOOD[ind1, ]
  c1_vect <- terra::vect(c1[, NAMES], geom = NAMES)
  er1 <- terra::rasterize(c1_vect, REF_RAST, fun = length)
  er1[is.na(er1)] <- 0
  err <- (er1 / er0)
  err[is.infinite(err)] <- 1
  err[err == 0] <- NA
  p <- terra::as.points(err)
  data <- cbind(terra::crds(p), terra::values(p))
  return(data)
}

gimmeAbsences <- function(IND0, YC_CROP, REF_RAST) {
  IND0 <- which(rowSums(YC_CROP[, sps]) == 0)
  cood0_vect <- terra::vect(YC_CROP[IND0, c("x", "y")], geom = c("x", "y"))
  er0 <- terra::rasterize(cood0_vect, REF_RAST, fun = length)
  p0 <- terra::as.points(er0)
  er0_df <- cbind(terra::crds(p0), terra::values(p0))
  er0_df$lyr.1 <- 0

  return(er0_df)
}

# Functions
makeSpCluster <- function(GENERA, DIST) {
  # I just want to make sure every species will be covered only once
  # Okay, we will go through the genera list alphabetically
  # Distance is set to 1 for all, except in the case of Androdon when distance
  # is 1.5 (obviously this is for a standardized tree)
  # We'll just make a list with the break up to begin with
  # Keep track of remaining genera

  # Initialize, list and counter
  genera2 <- GENERA
  spList <- list()
  df <- list()
  count <- 1
  nn <- ""

  for (i in seq_len(length(GENERA))) {
    genus <- GENERA[i]
    if (genus %in% genera2) {
      print(genus)
      ind <- grep(genus, rownames(DIST))

      # Set cut-offs
      if (genus == "Androdon") {
        cutoff <- 1.5
      } else {
        cutoff <- 1
      }
      sps <- spByDist(DIST, cutoff, ind[1])
      print(length(sps))
      spList[[count]] <- sps
      nn <- append(nn, genus)
      # List all selected genera
      inc <- getGenus(sps)
      leftBehindWarning(rownames(DIST), sps)

      # and remove from genera2
      rmv <- which(genera2 %in% inc)
      genera2 <- genera2[-rmv]
      count <- count + 1
    } else {
      next
    }
  }
  nn <- nn[-1]
  names(spList) <- nn
  return(spList)
}

getGenus <- function(SP) {
  genus <- substr(SP, 1, regexpr("_", SP) - 1)
  return(genus)
}

leftBehindWarning <- function(SP, SP_SEL) {
  SP2 <- SP[-which(SP %in% SP_SEL)]
  g2 <- getGenus(SP2)
  g <- getGenus(SP_SEL)
  w <- which(g %in% g2)
  if (length(w)) print(warning("Broken at genus level"))
}

spByDist <- function(DIST, CUTOFF, SPIND) {
  sel <- which(DIST[, SPIND] < CUTOFF)
  sps <- unique(c(names(sel), colnames(DIST[, SPIND])))
}

correctSpCluster <- function(SPLIST, GROUPS, INDEX) {
  # Set up
  count <- length(SPLIST) + 1

  # First loop over index where the spList needs to be broken down
  for (i in seq_len(length(INDEX))) {
    sps <- SPLIST[[INDEX[i]]]
    splits <- GROUPS[[i]]
    # Second loop over the length of the splits with each cluster
    for (j in seq_len(length(splits))) {
      if (j == 1) {
        SPLIST[[count]] <- sps[1:splits[j]]
      } else {
        SPLIST[[count]] <- sps[(splits[j - 1] + 1):splits[j]]
      }
      count <- count + 1
    }
  }
  return(SPLIST)
}

getEnvRasters <- function(ENV_CRS, EPATH, ENV_FILES, AREA_VECT = NULL, EX) {
  env_paths <- file.path(EPATH, ENV_FILES)
  env <- list()
  ann <- list()
  for (ii in seq_len(length(env_paths))) {
    env_path <- env_paths[ii]
    print(paste0("Loading ", env_path, " ..."))
    env_layer <- terra::rast(env_path)

    if (!is.null(AREA_VECT)) {
      env_mask <- terra::mask(env_layer, area_vect)
    }

    if (!is.null(EX)) {
      env_crop <- terra::crop(env_mask, ex)
    }
    env[[ii]] <- env_crop
  }
  return(env)
}


annotateCoods <- function(EPATH, ENV_FILES, ENV_CRS, COOD) {
  env_paths <- file.path(EPATH, ENV_FILES)
  cood <- COOD[, c("lon", "lat")]
  colnames(cood) <- c("x", "y")
  extent <- terra::ext(getExtent(cood))

  env <- list()
  for (i in seq_len(length(ENV_FILES))) {
    env_path <- env_paths[i]
    print(paste0("Loading ", env_path, " ..."))
    env_layer <- terra::rast(env_path)

    # Reproject
    # print(paste0(" Reprojecting..."))
    # env_layer <- terra::project(env_layer, ENV_CRS)

    if (!is.null(extent)) {
      env_crop <- terra::crop(env_layer, extent)
    }
    env[[i]] <- env_crop
  }
  cood_mat <- do.call(cbind, cood)

  ###### 6. Annotate coordinates with envs #######
  a <- list()
  for (i in seq_len(length(ENV_FILES))) {
    print(i)
    a[[i]] <- terra::extract(x = env[[i]], y = cood_mat)
  }

  annotations <- do.call(cbind, a)
  if (any(grepl("\\.tif", ENV_FILES))) {
    if (any(grepl("\\.tiff", ENV_FILES))) {
      xx <- gsub("\\.tiff", "", ENV_FILES)
    } else {
      xx <- gsub("\\.tif", "", ENV_FILES)
    }
  }
  colnames(annotations) <- xx
  return(annotations)
}


# EDA plots
makeEDA_pres_BG <- function(Y_TRAIN, X_TRAIN) {
  pdf(file = paste0("~/Downloads/test_", names(spList)[i], ".pdf"))
  par(mfrow = c(4, 2))
  for (j in seq_len(ncol(Y_TRAIN))) {
    for (k in seq_len(ncol(X_TRAIN))) {
      tryCatch(
        {
          abs_density <- max(density(X_TRAIN[which(Y_TRAIN[, j] == 0), k])$y)
          pres_density <- max(density(X_TRAIN[which(Y_TRAIN[, j] == 1), k])$y)
          ymax <- max(abs_density, pres_density)
          plot(density(X_TRAIN[which(Y_TRAIN[, j] == 0), k]),
            col = "red",
            main = paste(colnames(Y_TRAIN)[j], ": ", colnames(X_TRAIN)[k]),
            ylim = c(0, ymax)
          )
          lines(density(X_TRAIN[which(Y_TRAIN[, j] == 1), k]))
        },
        error = function(msg) {
          message(paste("Error for member ", j))
          # print(msg)
          return(NA)
        }
      )
    }
  }
  dev.off()
}

# rasterize cood
rasterizeCood <- function(YC_CROP, SPS, ER) {
  data <- list()
  for (j in seq_len(length(SPS))) {
    print(SPS[j])
    yfsp <- YC_CROP[, SPS[j]]
    # Get indices for presences/absences
    ind0 <- which(yfsp == 0)
    ind1 <- which(yfsp == 1)

    # Returns the ratio of presences to absences at every pixel
    err_df <- gimmeRatios(yfsp, YC_CROP[, c("x", "y")], ER)
    err_df$sp <- SPS[j]
    data[[j]] <- err_df
  }
  return(data)
}

makeTrainDataset <- function(Y, COOD, X, BG_COOD, Y_BG, X1) {
  # For each species, take out 30% of presences
  # Let absences be presences of other species
  # This can be modified though (by adding additional absences)

  # For the test set, each species should have that 30%
  # presences and an equivalent number of absences sampled
  # from neighboring cells
  # Loop to get presences indices for each species
  pres_ids <- list()
  for (k in seq_len(ncol(Y))) {
    pres_ids[[k]] <- which(Y[, k] == 1)
  }
  names(pres_ids) <- colnames(Y)

  pres_id_train <- list()
  pres_id_test <- list()
  for (k in seq_len(ncol(Y))) {
    n70 <- round(0.7 * length(pres_ids[[k]]))
    train_pres <- sample(1:length(pres_ids[[k]]), n70, replace = FALSE)
    pres_id_train[[k]] <- pres_ids[[k]][train_pres]
    pres_id_test[[k]] <- pres_ids[[k]][which(1:length(pres_ids[[k]])
    %in% train_pres == FALSE)]
  }
  # names(pres_id_train) <- names(pres_id_test) <- colnames(y)

  # Training dataset
  pres_id_train <- do.call(c, pres_id_train)
  pres_id_train <- pres_id_train[!duplicated(pres_id_train)]

  y_train <- Y[pres_id_train, ]
  cood_train <- COOD[pres_id_train, ]
  x_train <- X[pres_id_train, ]

  # Adding some background points
  # Take 70% of the background points for training and
  # 30% for testing
  bg_n_train <- 0.7 * nrow(BG_COOD)
  bg_train_ind <- sample(seq_len(nrow(BG_COOD)), bg_n_train)

  bg_cood_train <- BG_COOD[bg_train_ind, ]
  y_bg_train <- Y_BG[bg_train_ind, ]
  x_bg_train <- X1[bg_train_ind, ]

  # Rbind
  y_train <- rbind(y_train, y_bg_train)
  x_train <- rbind(x_train, x_bg_train)
  cood_train <- rbind(cood_train, bg_cood_train)

  # Test dataset
  cood_test_pool <- COOD[-pres_id_train, ]
  y_test_pool <- Y[-pres_id_train, ]
  x_test_pool <- X[-pres_id_train, ]

  # Now repeat background sampling for test dataset
  bg_cood_test <- BG_COOD[-bg_train_ind, ]
  y_bg_test <- Y_BG[-bg_train_ind, ]
  x_bg_test <- X1[-bg_train_ind, ]

  # Rbind
  y_test_pool <- rbind(y_test_pool, y_bg_test)
  x_test_pool <- rbind(x_test_pool, x_bg_test)
  cood_test_pool <- rbind(cood_test_pool, bg_cood_test)

  return(list(
    "y_train" = y_train, "x_train" = x_train,
    "cood_train" = cood_train, "y_test_pool" = y_test_pool,
    "x_test_pool" = x_test_pool, "cood_test_pool" = cood_test_pool,
    "pres_id_test" = pres_id_test, "pres_id_train" = pres_id_train,
    "pres_ids" = pres_ids
  ))
}


makeTestDataset <- function(
    COOD, EX, Y, X, PRES_ID_TEST, YC_TEST_POOL,
    PRES_IDS) {
  # Loop for test data per species to even out number of absences
  # per species in order not to skew eval metrics
  # I think test dataset should be per species?
  # It can be a list per species
  sp_test_data <- list()
  sp_test_cood <- list()
  sp_test_x <- list()

  for (k in seq_len(ncol(Y))) {
    # First get presences extent
    # Create modeling domain by increasing extent by X
    # Then add remaining 30% of presences to new dataframe
    # Sample equal number of absences from the domain extent
    # Add to presence dataframe
    # Add to a species list

    tmpy <- Y[PRES_ID_TEST[[k]], ]
    tmpcood <- COOD[PRES_ID_TEST[[k]], ]
    tmpx <- X[PRES_ID_TEST[[k]], ]

    yc_test_pool2 <- cropDFCood(YC_TEST_POOL, EX)
    rmv <- which(yc_test_pool2[, (k + 2)] > 0)
    if (length(rmv) > 0) yc_test_pool2 <- yc_test_pool2[-rmv, ]
    n <- length(PRES_ID_TEST[[k]])

    s0 <- seq_len(nrow(yc_test_pool2))
    if (n > length(s0)) n <- length(s0)
    ind0 <- sample(seq_len(nrow(yc_test_pool2)), n, replace = FALSE)
    tmpy0 <- yc_test_pool2[ind0, -c(1, 2)]
    tmpcood0 <- yc_test_pool2[ind0, c(2, 1)]
    colnames(tmpcood0) <- c("lat", "lon")
    tmpx0 <- x_test_pool[ind0, ]

    sp_test_data[[k]] <- rbind(tmpy, tmpy0)
    sp_test_cood[[k]] <- rbind(tmpcood, tmpcood0)
    sp_test_x[[k]] <- rbind(tmpx, tmpx0)
  }

  names(sp_test_data) <- names(sp_test_cood) <- colnames(Y)
  names(sp_test_x) <- colnames(Y)

  return(list(
    "sp_test_data" = sp_test_data, "sp_test_cood" = sp_test_cood,
    "sp_test_x" = sp_test_x
  ))
}

main <- function() {
  if (length(sps) == 1) next
  print(names(spList)[i])
  y1 <- raw[, sps]
  tree <- ape::keep.tip(tree1, tip = which(tree1$tip.label %in% sps))

  # This is the extent for this cluster
  pres_rw <- which(rowSums(y1) > 0)
  pres_cood <- cood1[pres_rw, ]
  ex <- getExtentDf(pres_cood, NAMES = c("lat", "lon"))

  # Widen extent
  ex <- terra::ext(ex + (c(-1, 1, -1, 1) * widenby))

  # Crop coords
  yc_crop <- cropDFCood(cbind(cood1, y1), ex)

  # This is the empty raster at 1km resolution
  er <- createSampleRast(ex, area_vect, epath, env_files[1])

  # Note to self - might not need the ratios method for the background
  # sampling scheme
  data <- rasterizeCood(yc_crop, sps, er) # uses the ratio method

  data_melt <- reshape2::melt(data, id = c("x", "y", "sp"))
  data_sps <- reshape2::dcast(data_melt, "x + y ~ sp", fun = sum)
  # For now, everything that is positive is awarded a 1
  data_sps[, sps] <- ifelse(data_sps[, sps] > 0, 1, 0)

  # Break into y and cood dataframes
  y <- data_sps[, sps]
  cood <- data_sps[, c(1, 2)]
  colnames(cood) <- c("lon", "lat")
  cood <- cood[, c(2, 1)]

  # Annotate coods
  x <- annotateCoods(epath, env_files, env_crs, COOD = cood)
  xx <- annotateCoods(
    "~/phylo-sdms/phyloproj/raw_data/hummingbirdSA",
    eff_files, env_crs, cood
  )
  colnames(xx) <- c("distance", "duration", "num_observers")
  x <- cbind(x, xx)

  # Complete cases only
  rmv <- which(complete.cases(x) == FALSE)
  if (length(rmv) > 1) {
    x <- x[-rmv, ]
    y <- y[-rmv, ]
    cood <- cood[-rmv, ]
  }

  # Remove all 0s
  rmv <- which(rowSums(y) == 0)
  if (length(rmv) > 1) {
    y <- y[-rmv, ]
    cood <- cood[-rmv, ]
    x <- x[-rmv, ]
  }
  # And TA-DA! That's the full dataset :)

  # Remove all rowsums in y that == 0
  # This is the baseline dataset where absences are
  # only presences of other species
  bg <- createBgPoints(cood, area_vect, ex, nrow(y), ncol(y))
  # Create vars
  y_bg <- bg$y_bg
  bg_cood <- bg$bg_cood
  x1 <- bg$x1

  # This right here (the rbind) might not be very useful
  colnames(y_bg) <- colnames(y)
  y <- rbind(y, y_bg)
  x <- rbind(x, x1)
  x[, "distance"] <- log(x[, "distance"])
  x[, "duration"] <- log(x[, "duration"])
  x[, "num_observers"] <- log(x[, "num_observers"])
  cood <- rbind(cood, bg_cood)

  # Saving this preemtively here
  cyx_sv <- cbind(cood, y, x)

  # makeOccurrenceMap(yc, FSP = sps[1])

  # Make the train and test set now
  train <- makeTrainDataset(y, cood, x, bg_cood, y_bg, x1)
  # Create vars
  y_train <- train$y_train
  x_train <- train$x_train
  cood_train <- train$cood_train
  x_test_pool <- train$x_test_pool
  cood_test_pool <- train$cood_test_pool
  pres_id_test <- train$pres_id_test
  pres_id_train <- train$pres_id_train
  pres_ids <- train$pres_ids

  # If you want to see stuff on a map---
  # yc_train <- cbind(cood_train, y_train)
  # makeOccurrenceMap(yc_train, sps[4])

  # Create test dataset
  # Test pool
  yc_test_pool <- cbind(cood_test_pool, y_test_pool)
  test <- makeTestDataset(
    cood, ex, y, x, pres_id_test,
    yc_test_pool, pres_ids
  )
  # Create vars
  sp_test_data <- test$sp_test_data
  sp_test_cood <- test$sp_test_cood
  sp_test_x <- test$sp_test_x

  # Write out all files
  store <- list(
    "y" = y_train,
    "x" = x_train,
    "cood" = cood_train
  )

  store_test <- list(
    "y" = sp_test_data,
    "x" = sp_test_x,
    "cood" = sp_test_cood
  )

  # Some plots
  makeEDA_pres_BG(y_train, x_train)

  # Filepaths to save
  filename_train <- paste0("FULL_ALL_run_files_rep1", ".Rdata")
  if (!dir.exists(file.path(dpath, names(spList)[i]))) {
    dir.create(file.path(dpath, names(spList)[i]))
  }
  save(store, file = file.path(dpath, names(spList)[i], filename_train))

  filename_test <- paste0("TEST-", filename_train)
  save(store_test, file = file.path(dpath, names(spList)[i], filename_test))
}

# Scale env vars in x and save means/sds for later prediction
scale_env_and_save <- function(x, env_vars, skip_vars = c("Intercept"), save_csv = NULL) {
  stopifnot(is.data.frame(x) || is.matrix(x))

  # Keep only columns that exist
  vars <- intersect(env_vars, colnames(x))
  if (!length(vars)) {
    message("No matching environmental variables found to scale.")
    return(list(x_scaled = x, scales = NULL))
  }

  # Exclude variables you never want to z-scale (e.g., Intercept)
  vars_to_scale <- setdiff(vars, skip_vars)
  if (!length(vars_to_scale)) {
    message("No variables to scale after excluding skip_vars.")
    scales_df <- data.frame(variable_name = vars, mean = NA_real_, sd = NA_real_)
    if (!is.null(save_csv)) utils::write.csv(scales_df, save_csv, row.names = FALSE)
    attr(x, "env_scale") <- scales_df
    return(list(x_scaled = x, scales = scales_df))
  }

  # Ensure numeric
  nn <- sapply(x[, vars_to_scale, drop = FALSE], is.numeric)
  if (!all(nn)) {
    stop(
      "Non-numeric columns in vars_to_scale: ",
      paste(vars_to_scale[!nn], collapse = ", ")
    )
  }

  # Compute means/sds (guarding against NA/zero SD)
  env_means <- vapply(x[, vars_to_scale, drop = FALSE], mean, numeric(1), na.rm = TRUE)
  env_sds <- vapply(x[, vars_to_scale, drop = FALSE], sd, numeric(1), na.rm = TRUE)
  bad_sd <- !is.finite(env_sds) | env_sds == 0
  if (any(bad_sd)) env_sds[bad_sd] <- 1

  # Apply scaling
  x[, vars_to_scale] <- sweep(x[, vars_to_scale, drop = FALSE], 2, env_means, "-")
  x[, vars_to_scale] <- sweep(x[, vars_to_scale, drop = FALSE], 2, env_sds, "/")

  # Build a tidy scales table for ALL env_vars (including skipped ones)
  scales_df <- data.frame(
    variable_name = vars,
    mean = NA_real_,
    sd = NA_real_,
    row.names = NULL
  )
  # fill for the ones we actually scaled
  scales_df$mean[match(names(env_means), scales_df$variable_name)] <- env_means
  scales_df$sd[match(names(env_sds), scales_df$variable_name)] <- env_sds
  # for skipped vars (e.g., Intercept) keep NA to signal "do not scale" later

  # Save to object attribute and optional CSV
  attr(x, "env_scale") <- scales_df
  if (!is.null(save_csv)) utils::write.csv(scales_df, save_csv, row.names = FALSE)

  message("Scaled env vars: ", paste(vars_to_scale, collapse = ", "))
  if (length(setdiff(vars, vars_to_scale))) {
    message("Skipped (not scaled): ", paste(setdiff(vars, vars_to_scale), collapse = ", "))
  }

  list(x_scaled = x, scales = scales_df)
}

# --- Use on prediction data later ---
# Applies the same means/sds; leaves vars with NA mean/sd untouched (e.g., Intercept).
apply_saved_scales <- function(x_pred, scales_df) {
  stopifnot(is.data.frame(x_pred) || is.matrix(x_pred))
  # Only vars that have finite mean/sd should be scaled
  scale_rows <- with(scales_df, is.finite(mean) & is.finite(sd) & sd != 0)
  vars_to_scale <- scales_df$variable_name[scale_rows]
  if (!length(vars_to_scale)) {
    return(x_pred)
  }

  # Ensure vars exist
  vars_to_scale <- intersect(vars_to_scale, colnames(x_pred))
  if (!length(vars_to_scale)) {
    return(x_pred)
  }

  means <- setNames(scales_df$mean[scale_rows], scales_df$variable_name[scale_rows])[vars_to_scale]
  sds <- setNames(scales_df$sd[scale_rows], scales_df$variable_name[scale_rows])[vars_to_scale]

  x_pred[, vars_to_scale] <- sweep(x_pred[, vars_to_scale, drop = FALSE], 2, means, "-")
  x_pred[, vars_to_scale] <- sweep(x_pred[, vars_to_scale, drop = FALSE], 2, sds, "/")
  x_pred
}
