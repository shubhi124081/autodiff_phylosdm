# condPred file
# Conditional prediction functions


#' (Compute X star)
#'
#' X star is used to calcuate covariance between focal sp and all other sp
#' The dimensions of X star is J-1, 1
#'
#' @param TREE_FULL full phylogenetic tree or distance matrix
#' @param J number of species
#' @param ALPHA mean of alpha posterior
#' @param RHO mean of rho posterior
#' @param FOCALSP focal species index number
#' @param SCALE scale phylogenetic tree/distance matrix or not (TRUE for yes)
#'
#' @returns a (J-1)x1 matrix that describes the relationship between the focal
#' species and all the other species
#'
#' @export

computeXstar <- function(TREE_FULL, SP, ALPHA, RHO, FOCALSP, SCALE = TRUE) {
    # Covar =  1 - t(Xstar) %*% X^-1 %*% Xstar

    L <- TREE_FULL
    if (SCALE) L <- L / max(L)
    x <- matrix(0, SP - 1, 1)

    for (i in seq_along(x)) {
        count <- i
        if (i >= FOCALSP) count <- i + 1
        x[i, 1] <- ALPHA^2 * exp(-0.5 * (L[count, 1])^2 / RHO^2)
    }

    return(x)
}


#' (Compute Xinv)
#'
#' X inv is the inverse of the covariance matrix of all species apart from the
#' focal species
#'
#' @param TREE_EXP phylogenetic tree or distance matrix used in the experiment
#' (tree_exp)
#' @param SP number of species
#' @param ALPHA mean of alpha posterior
#' @param RHO mean of rho posterior
#' @param SCALE scale phylogenetic tree/distance matrix or not (TRUE for yes)
#' @param ST_DEVS standard deviation vector
#'
#' @returns a JxJ matrix that describes the inverse covariance matrix of
#' species in the analysis
#'
#' @export


computeXinv <- function(TREE_EXP, SP, ALPHA, RHO, SCALE = TRUE, ST_DEVS) {
    L <- TREE_EXP
    if (SCALE) L <- L / max(L)

    x <- matrix(0, SP, SP)

    for (i in 1:(SP - 1)) {
        x[i, i] <- 1 + 0.1 # for numerical stability
        for (j in (i + 1):SP) {
            x[i, j] <- ALPHA^2 * exp(-0.5 * (L[i, j])^2 / RHO)
            x[j, i] <- x[i, j]
        }
    }

    for (n in seq_along(SP)) x[n, n] <- 1 + 0.1
    diag(x) <- diag(x) * ST_DEVS
    xinv <- solve(x)
    return(x)
}

#' check and (Convert Phylo) object to a distance (Matrix)
#'
#' Check whether tree has been loaded in as a phylogenetic tree (i.e. object of
#' class phylo) and if so, it needs to be converted to a co-phenetic distance
#' matrix object. If it is already loaded as a co-phenetic distance matrix,
#' nothing is done
#'
#' @param TREE tree object to check and convert
#'
#' @returns a co-phenetic distance matrix of dimensions JxJ
#'
#' @export

convertPhyloMatrix <- function(TREE) {
    if (class(TREE) == "phylo") {
        TREE <- ape::cophenetic.phylo(TREE)
        print("phylo object converted to distance matrix")
    } else {
        print("already a distance matrix")
    }

    return(TREE)
}


#' Get model posteriors
#'
#' Grab model posteriors given user specified variables. Option to return the
#' whole posterior or aggregate by a function for example, mean or standard dev
#'
#' @param VAR variable to return
#' @param N n is a range of vectors e.g. 1 to K
#' @param SSPEC species specific or not
#' @param CSPEC covariate specific or not
#'
#' @return a list of posterior or posterior means/sds
#'
#' @export


getPosteriors <- function(VAR = "beta", N = c(1:K),
                          SSPEC = TRUE, CSPEC = TRUE) {
    var <- list()

    if (SSPEC && CSPEC) {
        for (i in N) {
            p <- paste0(VAR, "[", i, ",", 1:J, "]")
            ind <- which(colnames(fit) %in% p)
            var[[i]] <- fit[, ind]
        }
    } else {
        if (SSPEC || CSPEC) {
            p <- paste0("^", VAR, "\\[")
            ind <- grep(p, colnames(fit))
            var <- fit[, ind]
        } else {
            ind <- which(colnames(fit) %in% VAR)
            var <- fit[, ind]
        }
    }

    return(var)
}

#' Compute conditional predictions
#'
#' @param FSPALL (F)ocal (S)pecies (ALL) all the species to do conditional
#' prediction on
#' @param TREE_FULL the entire phylogenetic distance matrix JXJ
#' @param SP number of species in total
#' @param TREE_EXP_FLAG the experimental tree - without the missing species.
#' If NULL, computed iteratively (J-1)x(J-1)
#' @param RHO rho mean posterior
#' @param ALPHA alpha mean posterior
#' @param NITER number of iterations
#' @param VARS total number of predictor variables
#' @param VFLAG variable flag - if "ENV" these are environmental predictors
#' @param WOPHYLO phylogenetic conditional prediction or not (W)ith(o)ut
#' (phylo)? If TRUE, independent cond pred
#' @param BETA list of beta posterior means
#' @param ST_DEVS vector of standard deviations for betas
#' @param WRITEOUT write out csv samples or not
#' @param WODIR write out directory
#' @param PRINT print out statements (TRUE) or not (FALSE)
#' @param DATADEF Is this a data deficient experiment or not
#' (yes/TRUE, no/FALSE)
#' @param EXP_FLAG which experiment are you running e.g. ALL, NOCC, NOPHYLO,
#' DEF5, etc
#' @param GENUS the focal genus
#' @param REP the rep number
#' @param DEFLEV the deficiency level [nodef, y100, y25 etc]. Except for nodef
#' the deflevel should be followed by y and the function assumes in the file
#' name, the deficiency level is surrounded by punctuation e.g. _y50_
#'
#' @return Nothing returned, files are written out
#'
#' @export


computeCondPred <- function(FSPALL = 1:J,
                            TREE_FULL = tree_full,
                            SP = J,
                            TREE_EXP = NULL,
                            RHO = rho,
                            ALPHA = alpha,
                            NITER = 3500,
                            VARS = 1:K,
                            VFLAG = "ENV",
                            WOPHYLO,
                            BETA = beta,
                            ST_DEVS = st_devs,
                            WRITEOUT = TRUE,
                            WODIR = paste0(res.directory, "/condPred_res"),
                            PRINT = TRUE,
                            DATADEF = FALSE,
                            EXP_FLAG,
                            GENUS = genus,
                            REP = rep,
                            DEFLEV = deflev) {
    # if(WEFFORT) {
    #   V <- (KONLY+1):VARS
    # } else {
    #   V <- 1:KONLY
    # }

    # Set up first loop - iterates over focal species
    sds <- matrix(0, nrow = length(FSPALL), ncol = 1)
    nm <- colnames(TREE_FULL)[FSPALL]
    pb <- txtProgressBar()

    for (o in seq_along(FSPALL)) {
        focalsp <- FSPALL[o]

        # If TREE_EXP is null, the experiment is to be run iteratively
        if (is.null(TREE_EXP)) {
            TREE_EXP <- TREE_FULL[-focalsp, -focalsp]
            st_devs <- ST_DEVS[-focalsp]
        }

        # Compute components of the covariance function
        if (WOPHYLO) {
            X <- diag(1, nrow = SP - 1, ncol = SP - 1)
            if (VFLAG == "ENV") diag(X) <- st_devs
            Xinv <- solve(X)
            Xstar <- diag(1, nrow = SP - 1, ncol = 1)
        } else {
            Xstar <- computeXstar(TREE_FULL, SP, ALPHA, RHO, focalsp)
            Xinv <- computeXinv(TREE_EXP,
                SP = nrow(TREE_EXP), ALPHA, RHO,
                ST_DEVS = st_devs
            )
        }

        covarFun <- t(Xstar) %*% Xinv %*% Xstar
        # should it be 1 - t(Xstar)%*%Xinv %*% Xstar
        sdEst <- sqrt(covarFun)
        sds[o, ] <- sdEst

        # Set up second loop - iterates over variables of interest
        means <- list()
        samples <- matrix(0, NITER, length(VARS))

        for (i in seq_along(VARS)) {
            varno <- VARS[i]

            b <- BETA[[varno]]
            if (length(b) == nrow(TREE_FULL)) {
                fspb <- b[which(names(b) == paste0("beta[", varno, ",", focalsp, "]"))]
                ospb <- b[which(names(b) != paste0("beta[", varno, ",", focalsp, "]"))]
            } else {
                ospb <- b
            }

            if (DATADEF) {
                means[[i]] <- 0 + t(Xstar) %*% Xinv %*% ospb # Removed fspb from sum
            } else {
                means[[i]] <- t(Xstar) %*% Xinv %*% ospb
            }

            samples[, i] <- rnorm(NITER, means[[i]], sdEst)
        }

        # Print info or progress bar
        if (PRINT) {
            print(paste0("Species no: ", o))
            print(colMeans(samples))
            bb <- paste0("beta[", VARS, ",", focalsp, "]")
            if (length(bb) == 1) {
                print(mean(fit[, bb]))
            } else {
                print(colMeans(fit[, bb]))
            }
        } else {
            setTxtProgressBar(pb, value = (o / length(FSPALL)))
        }

        # Write out files - change name depending on effort/or not
        if (WRITEOUT) {
            np <- file.path(WODIR, GENUS)
            dir.create(np)
            colnames(samples) <- paste0("beta_", VARS)


            if (WOPHYLO) {
                write.csv(samples,
                    file = file.path(np, paste0(
                        EXP_FLAG, "_", nm[o],
                        "_", "cond_samples_", VFLAG, "_", DEFLEV, "_", "rep", REP, ".csv"
                    )),
                    row.names = FALSE
                )
                # write.csv(sds, file =file.path(np, paste0(EXP_FLAG, nm[o],"_",
                # "sd_est_EFFORT.csv")), row.names = F)
            } else {
                write.csv(samples,
                    file = file.path(np, paste0(
                        EXP_FLAG, "_", nm[o],
                        "_", "cond_samples_", VFLAG, "_", DEFLEV, "_", "rep", REP, ".csv"
                    )),
                    row.names = FALSE
                )
                # write.csv(sds, file =file.path(np, paste0(EXP_FLAG, nm[o],"_",
                # "sd_est_ENV.csv")), row.names = F)
            }
        }
    }
}


#' (Pick) (F)ocal (SP)ecies (Beta)
#'
#' Reformat beta posteriors s.t. a list of beta posteriors is converted
#' to a matrix of betas for a focal species
#'
#' @param BETA list of beta posteriors with length = number of betas and each
#' element NITER x J where NITER is number of posterior iterations and J is
#' number of species
#' @param FSPIND Focal species index - a number
#'
#' @return A matrix of beta posteriors for the focal species
#'
#' @export

pickFSPBeta <- function(BETA, FSPIND) {
    xx <- do.call(cbind, BETA)

    i <- grep(paste0(",", FSPIND, "\\]"), colnames(xx))
    w3 <- xx[, i]

    return(w3)
}


#' Little worker function to get the deficiency level from an experiment's name
#'
#' @description there are a couple of rules to take note of here, first of all,
#' the experiment name needs to have deficiency level in the name and it needs
#' to be wrapped by some sort of punctuation, most likely underscores. Second,
#' if the experiment is not data-deficient, "nodef" is the tag that needs to
#' be included in the experiment name
#'
#' @param EXP_NAME experiment name
#'
#' @return a string describing the deficiency level or an error in case the
#' function is unable to find the deficiency level
#'
#' @export

getDefLev <- function(EXP_NAME) {
    # Get deficiency level
    deflev <- ""
    if (grepl("nodef", EXP_NAME)) {
        deflev <- "nodef"
    } else {
        tmp <- regmatches(EXP_NAME, gregexpr("[[:punct:]]y[0-9]{1,3}", EXP_NAME,
            perl = TRUE
        ))
        deflev <- substr(tmp, 2, nchar(tmp))
    }
}


# Plotting functions ------

# Functions for plotting figures
# dummy edit
#' Converting a species name into a species code
#' A species code is defined as the first two letters of the genus and first
#' two letters of the species name. This can be modified if there are duplicates
#' (sometimes there are)
#'
#' @param X the species name - the genus and species have to be separated
#' by a "_"
#' @param CUTOFF the cutoff for species name - 2 letters, 3, letters etc.
#' Often 4 letters ensures no duplicates
#'
#' @export

toSpCode <- function(X, CUTOFF) {
    code1 <- substr(X, 1, 2)
    code2 <- substr(X, regexpr("_", X) + 1, regexpr("_", X) + CUTOFF)
    code2 <- paste(toupper(substr(code2, 1, 1)),
        substr(code2, 2, nchar(code2)),
        sep = ""
    )
    shrtname <- paste(code1, code2, sep = "")
    return(shrtname)
}




# A phylogenetic tree plot with an accompanying bar plot to show number of
# presence records per species
#
# Creates a pot of the phylogenetic tree and a bar plot of number of
# presence records
# NOTE: toSpCode needs to be packaged in a general package and called
# For now, it needs to be in the general environment, sourced from
# ~/Documents/R/shubhiFunctions.R
#
# @param TREE the tree (this is a phylo object)
# @param Y the occurrence records data frame for the species in tree
# @param LAYOUT layout of the tree - common options are rectangular or circular.
# For more see ggtree documentation
#
# @importFrom ggtreeExtra geom_fruit
# @import ggplot2
# @import ggtree
#
# @return A ggplot object (map in this case)
#
# @export


plotOccurrenceTree <- function(TREE, Y, LAYOUT = "rectangular", CUTOFF = 2) {
    library(ggtree)
    library(ggplot2)

    shrtname <- toSpCode(TREE$tip.label, CUTOFF)
    ww <- which(duplicated(shrtname))
    shrtname[ww] <- paste0(shrtname[ww], "2")

    TREE$tip.label <- shrtname

    trplot <- ggtree(TREE, layout = LAYOUT) + geom_tiplab(size = 2)

    df <- data.frame("total" = colSums(Y))
    rwns <- toSpCode(colnames(Y), CUTOFF)
    ww <- which(duplicated(rwns))
    rwns[ww] <- paste0(rwns[ww], "2")
    rownames(df) <- rwns

    codnms <- toSpCode(colnames(Y), CUTOFF)
    ww <- which(duplicated(codnms))
    codnms[ww] <- paste0(codnms[ww], "2")
    df$code <- codnms

    df$category <- ifelse(df$total > 30, "Rich", "Deficient")
    trbar <- trplot + ggtreeExtra::geom_fruit(
        data = df,
        geom = geom_bar,
        mapping = aes(
            x = total,
            y = code,
            fill = category
        ),
        orientation = "y",
        stat = "identity",
        pwidth = 1,
        offset = .2
    ) +
        scale_fill_manual(values = c("grey", "dodgerblue4")) +
        labs(fill = "") +
        ggtitle("Phylogenetic tree and number of presence records",
            subtitle = paste0("Note: Species with fewer than 30 presence
points are categorized as data-deficient
                              \nMaximum number of records =
", max(df$total), "\nMinimum number of records =", min(df$total))
        )

    return(trbar)
}


#' Get Environmental variables for a specified extent
#' @param RASTER2 raster stack data to be returned
#' @param COOD2 coordinates (used to calculate extent).
#' Must be labelled "x", "y" test test
#'
#' @return cropped raster stack
#'
#' @importFrom raster crop
#' @importFrom raster xyFromCell
#' @importFrom raster ncell
#' @importFrom raster getValues
#'
#' @export

# getEnv <- function(COOD2, RASTER2) {
#     ex <- terra::ext(COOD2)
#     xmin <- ex[1]
#     xmax <- ex[2]
#     ymin <- ex[3]
#     ymax <- ex[4]

#     ras <- terra::crop(RASTER2, ex)


#     xcoords <- raster::xyFromCell(ras, seq_len(raster::ncell(ras)))
#     vals <- raster::getValues(ras)
#     rasdf <- data.frame(cbind(xcoords, vals))
#     colnames(rasdf) <- c("x", "y", "value")
#     rasdf <- rasdf[-c(which(is.na(rasdf$value))), ]

#     return(rasdf)
# }

#' Make Environment variable Map
#'
#' Creates a plot of the given enviornmental variable. To overlay a data frame,
#' for example,
#' a presence point data frame
#' NOTE: scale_fill_colorway_c needs to be packaged in RColorways and loaded in
#' as a library,
#' For now, it needs to be in the general environment, sourced from
#' ~/color_pkg/main.R
#'
#' @param RASTER raster to be plotted
#' @param COOD coordinates (used to calculate extent). Must be labelled "x", "y"
#' @param RTITLE Raster title to be printed on map
#' @param YC (optional) The occurrence record
#' @param FSP (F)ocal (sp)ecies name
#' @param PALETTE color palette (a hex code vector) to be used to color the map
#'
#' @importFrom rnaturalearth ne_countries
#' @importFrom raster extent
#' @importFrom raster xyFromCell
#' @importFrom raster crop
#' @importFrom raster getValues
#' @importFrom raster ncell
#'
#' @import ggplot2
#'
#' @return A ggplot object (map in this case)
#'
#' @export

makeEnvMap <- function(RASTER, COOD, RTITLE, YC = NULL, FSP = NULL,
                       PALETTE = ALL$Divergent$`Rd-Y-Bl`, WORLD = world) {
    warning("This plot takes 2-3 minutes")


    rasdf <- getEnv(COOD, RASTER)
    ex <- terra::ext(COOD)
    xmin <- ex[1]
    xmax <- ex[2]
    ymin <- ex[3]
    ymax <- ex[4]

    if (is.null(YC)) {
        plt <- ggplot(rasdf) +
            geom_tile(aes(x = x, y = y, fill = value)) +
            geom_sf(data = WORLD, fill = NA, col = "black") +
            coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) +
            theme_bw() +
            theme(panel.grid = element_blank()) +
            scale_fill_colorway_c(palette = PALETTE) +
            ggtitle(paste0("Map of ", RTITLE))
    } else {
        # Overlay
        yfsp <- YC[, FSP]
        ind <- which(yfsp == 1)
        yc2 <- YC[ind, c(1, 2)]

        plt <- ggplot(rasdf) +
            geom_tile(aes(x = x, y = y, fill = value)) +
            geom_sf(data = world, fill = NA, col = "black") +
            geom_point(
                data = yc2, aes(x = x, y = y), color = "black",
                shape = 3
            ) +
            coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) +
            theme_bw() +
            theme(panel.grid = element_blank()) +
            scale_fill_colorway_c(palette = PALETTE) +
            ggtitle(paste0("Map of ", RTITLE),
                subtitle = paste0("Overlaid presence points for species ", FSP)
            )
    }

    return(plt)
}

#' Annotate presence/absence points with given a selected variable
#'
#' @param YXC Y data, X data and coods
#' @param VAR variable to annotate with e.g. temperature
#' @param PALETTE selected color scheme
#'
#' @return return annotated plot
#'
#' @import ggplot2
#'
#' @export


annotatePointsPlot <- function(YXC, VAR, FSP,
                               PALETTE = ALL$Divergent$`Rd-Y-Gn`) {
    varno <- which(colnames(yxc) == VAR)
    fspno <- which(colnames(yxc) == FSP)
    p <- ggplot(YXC) +
        geom_point(aes(
            x = YXC[, varno],
            y = as.factor(YXC[, fspno]),
            fill = YXC[, varno]
        ), shape = 21, size = 3) +
        scale_fill_colorway_c(palette = PALETTE) +
        theme_bw() +
        ylab("Presence/Absence") +
        xlab(VAR) +
        labs(fill = VAR) +
        ggtitle(paste0("Annotated plot for ", FSP),
            subtitle = paste0("Annotated with ", VAR)
        )
    return(p)
}



# (Paint) the (branches) of a tree with number of (occ)urrence records
#
# @param VALUES names vector of occurrence records with species names
# @param TREE the phylogenetic tree
# @param COLORS color palette
# @param RETURNDF return data frame
#
# @return a list with a plot with branches colored with category each
# species is depending on the number of presence records it has and optionally
# a dataframe put together within the function
#
# @import ggtree
# @importFrom ggtree ggtree
# @importFrom ggtree groupOTU

# @export


paintBranchesOCC <- function(VALUES, TREE, COLORS, RETURNDF = TRUE) {
    b <- quantile(VALUES, seq(0, 1, by = 0.1))
    b[1] <- 0

    df <- data.frame(
        "sp" = names(VALUES),
        "occ" = VALUES
    )

    df$int <- cut(VALUES, b)
    dict <- data.frame("int" = unique(df$int)[order(unique(df$int))])
    dict$cat <- paste0(seq(10, 100, by = 10), "%")
    dict$cat[nrow(dict)] <- "99%"

    df2 <- merge(df, dict, by = "int")
    rownames(df2) <- df2$sp

    m <- unique(df2$cat)
    m <- m[order(m)]
    grp <- list()
    for (i in 1:length(m)) grp[[i]] <- df2$sp[df2$cat == m[i]]

    names(grp) <- m

    # Occurrence painted on tree
    p <- ggtree::ggtree(TREE, layout = "circular")
    p1 <- ggtree::groupOTU(p, grp, "Occurrence") + aes(color = Occurrence) +
        theme(legend.position = "right") +
        scale_color_manual(values = cols)

    if (RETURNDF) {
        o <- list("plot" = p1, "df" = df2)
    } else {
        o <- list("plot" = p1)
    }

    return(o)
}


#' Read in multiple csv files
#'
#' Read in multiple csv files given names of selected files and
#' the directory
#'
#' @param DPATH (path) to (data)
#' @param files name of files - can have ".csv" in them
#'
#' @return read in data, named by file
#'
#' @export

readMultiple <- function(DPATH, FILES) {
    df <- data.frame("no" = seq(1, length(FILES)), "files" = FILES)
    store <- list()

    store <- apply(df, 1, function(x) {
        w <- read.csv(file.path(DPATH, x["files"]))
        return(w)
    })

    if (grepl(".csv", df$files[1])) {
        names(store) <- gsub(".csv", "", df$files)
    } else {
        names(store) <- df$files
    }

    return(store)
}

#' Simple function to write out the codes from the experiment name
#'
#' The experiment name variable must be "Var2"
#'
#' @param DF the dataframe with Var2 as the experiment column
#'
#' @return Added codes to DF
#'
#' @export


makeCodes <- function(DF) {
    tmp_code <- apply(DF, 1, function(x) {
        ww <- x["Var2"]
        getDefLev(ww)
    })
    DF$code <- tmp_code
    DF$code <- ifelse(grepl("nodef", DF$code), "FULL", DF$code)
    # DF$code <- ifelse(grepl("y25+|y5+", DF$Var2), DF$code, "FULL")
    # DF$code <- ifelse(grepl("y5", DF$Var2), "y5", DF$code)

    DF$exp <- ifelse(grepl("BASE", DF$Var2), "OCC", "PHYLO")
    DF$exp <- ifelse(grepl("_ALL", DF$Var2), "OCC+PHYLO", DF$exp)
    DF$exp <- ifelse(grepl("MISS", DF$Var2), "PHYLO", DF$exp)
    DF$exp <- ifelse(grepl("NOPHYLO", DF$Var2), "OCC", DF$exp)


    exps <- unique(DF$exp)
    if ("NO PHYLO" %in% exps) {
        warning("NO PHYLO shouldn't be in the experiments list")
    }

    DF$frank <- paste0(DF$code, "-", DF$exp)

    return(DF)
}



#' Insert experiment codes for AUC plot
#'
#' This is a nuisance function, jsut for putting in codes for the boxplot
#' Usually there are codes for level of data-deficiency (FULL, y25, y5)
#' and experiment (NOPHYLO, ALL, MISS). Both BASE and NOPHYLO count as NOPHYLO
#' since neither uses phylo information to conditionally predict parameters
#'
#' @param STORE a list or other that has all the data for the plot
#'
#' @return df to be used for the AUC plot
#'
#' @export

AUCcodes <- function(STORE) {
    # df <- reshape2::melt(STORE)
    expno <- length(unique(STORE$Var2))
    reps <- nrow(STORE) / expno

    df <- makeCodes(DF = STORE)

    # Create a dummy space
    codes <- paste0(unique(df$code), "_DUM")
    dummy <- df[1:(reps * length(codes)), ]
    dummy[, c(1, 2, 3)] <- rep(NA, nrow(dummy))
    dummy$code <- rep(codes, each = reps)
    df <- rbind(df, dummy)
    df$frank <- paste0(df$code, "-", df$exp)

    # dummy <- df[1:300, ]
    # dummy[, c(1, 2, 3)] <- rep(NA, nrow(dummy))
    # dummy$code <- rep(c("FULL", "y25", "y5"), each = reps)
    # dummy$exp <- rep("x", nrow(dummy))
    # df <- rbind(df, dummy)
    # df$frank <- paste0(df$code, "-", df$exp)

    return(df)
}


#' add jitter to the boxplot
#'
#' Another nuisance function, add a little nugget to the boxplot values to
#' make them more obvious on the plot
#'
#' @param DF the dataframe for the boxplot
#' @param IND the indices for experiments to jitter
#' @param NREPS number of AUC values
#'
#' @return return jittered DF
#'
#' @export

addJitter <- function(DF, IND, NREPS = 100) {
    for (i in seq_along(IND)) {
        ii <- IND[i]

        if (ii == 1) {
            DF[1:NREPS, "value"] <- DF[1:NREPS, "value"] + rnorm(NREPS, 0, 0.01)
        } else {
            range <- c(((ii * NREPS) + 1):((ii + 1) * NREPS))
            DF[range, "value"] <- DF[range, "value"] + rnorm(NREPS, 0, 0.01)
        }
    }

    return(DF)
}

#' Create Performance boxplot
#'
#' Create performance boxplots with specs for x, y etc for the ggplot
#'
#' @param DF data frame
#' @param Y y axis
#' @param X x axis
#' @param FILL what to fill by
#' @param LAB fill title
#' @param COLS color scheme
#' @param FSP focal performance
#' @param RANGE ylim range
#'
#' @return ggplot boxplot
#'
#' @export

performanceBoxPlot <- function(DF, Y, X, FILL, LAB, COLS, FSP, RANGE) {
    p <- ggplot(DF) +
        geom_hline(yintercept = 0.5, color = "grey60", linetype = "dashed") +
        geom_boxplot(aes(y = Y, x = X, fill = FILL)) +
        theme_bw() +
        theme(
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            panel.grid = element_blank()
        ) +
        ylim(RANGE) +
        ylab(LAB) +
        scale_fill_manual(values = COLS) +
        labs(fill = "") +
        ggtitle(paste0("AUC plots for ", FSP),
            subtitle = "Comparison of model performance"
        )

    return(p)
}


#' compute an improvement in AUC values compared to a standard
#'
#' @param COMP compare df has to have codes FULL, y25, y5
#' @param STANDARD standard df has to have codes FULL, y25, y5
#'
#' @return return new comp
#'
#' @export

computeImprovement <- function(COMP, STANDARD) {
    COMP$improvement <- rep(0, nrow(COMP))
    COMP$improvement <- c(
        COMP$value[COMP$code == "FULL"] -
            STANDARD$value[STANDARD$code == "FULL"],
        COMP$value[COMP$code == "y25"] -
            STANDARD$value[STANDARD$code == "y25"],
        COMP$value[COMP$code == "y5"] -
            STANDARD$value[STANDARD$code == "y5"]
    )

    return(COMP)
}

#' A sanity check of centroids, grids and range
#'
#' @param GRID Regional grid (terra vect)
#' @param RANGE Species range (terra vect)
#' @param CENTROIDS grid centroids (terra vect)
#' @param CENTROIDS2 range centroids (terra vect)
#'
#' @return a plot
#'
#' @importFrom terra plot
#' @importFrom terra lines
#' @importFrom terra points
#'
#' @export

sanityCheckGrid <- function(GRID, RANGE, CENTROIDS, CENTROIDS2) {
    terra::plot(GRID)
    terra::lines(RANGE)
    terra::points(CENTROIDS, col = "red")
    terra::points(CENTROIDS2, col = "blue")
}

#' Presence probability map
#'
#' @param DF dataframe containing coordinates and presence probability
#' (respectively named x, y, and pp)
#' @param COLS color scheme
#' @param TITLE title
#' @param OCC (logical) plot observed data or not
#' @param OCC_DF dataframe containing observed data (with x and y columns)
#' @param RANGE (logical) plot species expert range or not
#' @param RANGE_VECT species expert range as a terra vect object
#' @param EX extent vector with xmin, xmax, ymin and ymax
#'
#' @return a map
#'
#' @import ggplot2
#'
#' @export


presenceMap <- function(DF, COLS, VAR, TITLE = "", OCC, OCC_DF = NULL, RANGE,
                        RANGE_VECT = NULL, EX) {
    xmin <- EX[1]
    xmax <- EX[2]
    ymin <- EX[3]
    ymax <- EX[4]

    ind <- which(colnames(DF) == VAR)
    maxp <- max(DF[, ind])

    # Plot of only presence probability
    p <- ggplot(data = DF) +
        geom_tile(aes(x = x, y = y, fill = DF[, ind])) +
        geom_sf(data = world, fill = NA, col = "black") +
        scale_fill_gradientn(colors = COLS, limits = c(0, maxp)) +
        coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) +
        theme_bw() +
        theme(panel.grid = element_blank()) +
        labs(fill = "Presence Probability") +
        ggtitle(TITLE)

    # Plot of presence probability and species range
    if (RANGE) {
        p <- ggplot(data = DF) +
            geom_tile(aes(x = x, y = y, fill = DF[, ind])) +
            geom_sf(data = RANGE_VECT, fill = NA, col = "green") +
            geom_sf(data = world, fill = NA, col = "black") +
            scale_fill_gradientn(colors = COLS, limits = c(0, maxp)) +
            coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) +
            theme_bw() +
            theme(panel.grid = element_blank()) +
            labs(fill = "Presence Probability") +
            ggtitle(TITLE)
    }

    # Plot of presence probability and observed data
    if (OCC) {
        p <- ggplot(data = DF) +
            geom_tile(aes(x = x, y = y, fill = DF[, ind])) +
            geom_sf(data = world, fill = NA, col = "black") +
            scale_fill_gradientn(colors = COLS, limits = c(0, maxp)) +
            geom_point(data = OCC_DF, aes(x = x, y = y), size = .5) +
            coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) +
            theme_bw() +
            theme(panel.grid = element_blank()) +
            labs(fill = "Presence Probability") +
            ggtitle(TITLE)
    }

    # Plot of presence probability, observed data and species range
    if (RANGE && OCC) {
        p <- ggplot(data = DF) +
            geom_tile(aes(x = x, y = y, fill = DF[, ind])) +
            geom_sf(data = RANGE_VECT, fill = NA, col = "green") +
            geom_sf(data = world, fill = NA, col = "black") +
            scale_fill_gradientn(colors = COLS, limits = c(0, maxp)) +
            geom_point(data = OCC_DF, aes(x = x, y = y), size = .5) +
            coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) +
            theme_bw() +
            theme(panel.grid = element_blank()) +
            labs(fill = "Presence Probability") +
            ggtitle(TITLE)
    }

    return(p)
}


#' Presence binary probability map
#'
#' @param DF dataframe containing coordinates where binary predictions are 1
#' (respectively named x, y, and pp)
#' @param COLS color scheme
#' @param TITLE title
#' @param OCC (logical) plot observed data or not
#' @param OCC_DF dataframe containing observed data (with x and y columns)
#' @param RANGE (logical) plot species expert range or not
#' @param RANGE_VECT species expert range as a terra vect object
#' @param EX extent vector with xmin, xmax, ymin and ymax
#'
#' @return a map
#'
#' @import ggplot2
#'
#' @export


presenceMapBinary <- function(DF, COLS, VAR, TITLE = "", OCC, OCC_DF = NULL,
                              RANGE, RANGE_VECT = NULL, EX) {
    xmin <- EX[1]
    xmax <- EX[2]
    ymin <- EX[3]
    ymax <- EX[4]

    # ind <- which(colnames(DF) == VAR)
    # maxp <- max(DF[, ind])

    # Plot of only presence probability
    p <- ggplot(data = DF) +
        geom_point(aes(x = x, y = y), col = "red", size = .1) +
        geom_sf(data = world, fill = NA, col = "black") +
        # scale_fill_gradientn(colors = COLS, limits = c(0, maxp)) +
        coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) +
        theme_bw() +
        theme(panel.grid = element_blank()) +
        labs(fill = "Presence Probability") +
        ggtitle(TITLE)

    # Plot of presence probability and species range
    if (RANGE) {
        p <- ggplot(data = DF) +
            geom_point(aes(x = x, y = y), col = "red", size = .1) +
            geom_sf(data = RANGE_VECT, fill = NA, col = "green") +
            geom_sf(data = world, fill = NA, col = "black") +
            # scale_fill_gradientn(colors = COLS, limits = c(0, maxp)) +
            coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) +
            theme_bw() +
            theme(panel.grid = element_blank()) +
            labs(fill = "Presence Probability") +
            ggtitle(TITLE)
    }

    # Plot of presence probability and observed data
    if (OCC) {
        p <- ggplot(data = DF) +
            geom_point(aes(x = x, y = y), col = "red", size = .1) +
            geom_sf(data = world, fill = NA, col = "black") +
            # scale_fill_gradientn(colors = COLS, limits = c(0, maxp)) +
            geom_point(data = OCC_DF, aes(x = x, y = y), size = .5) +
            coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) +
            theme_bw() +
            theme(panel.grid = element_blank()) +
            labs(fill = "Presence Probability") +
            ggtitle(TITLE)
    }

    # Plot of presence probability, observed data and species range
    if (RANGE && OCC) {
        p <- ggplot(data = DF) +
            geom_point(aes(x = x, y = y), col = "red", size = .1) +
            geom_sf(data = RANGE_VECT, fill = NA, col = "green") +
            geom_sf(data = world, fill = NA, col = "black") +
            # scale_fill_gradientn(colors = COLS, limits = c(0, maxp)) +
            geom_point(data = OCC_DF, aes(x = x, y = y), size = .5) +
            coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) +
            theme_bw() +
            theme(panel.grid = element_blank()) +
            labs(fill = "Presence Probability") +
            ggtitle(TITLE)
    }

    return(p)
}


# postProcessing file -----


#' Extract environmental variables within a range
#' @description For a given expert range (or more generically, a terra vect)
#' file and a list of environmental variables, returns either a cropped raster
#' or a dataframe of extracted values within the range.
#' Note: epath is the filepath to the environmental variables, this needs to
#' be in the environment before the function is run
#'
#' @param ENV_FILE list with names of the environmental variables (tiff files)
#' @param FACT factor to aggregate by (default is 10)
#' @param RANGE species expert range as a terra vect file
#' @param PTS (logical) if TRUE return values else return raster file
#'
#' @return Depending on the value of PTS, either a dataframe or a raster
#'
#' @importFrom terra rast
#' @importFrom terra aggregate
#' @importFrom terra crop
#' @importFrom terra extract
#'
#' @export

extractEnvInRange <- function(ENV_FILE, FACT = 10, RANGE, PTS) {
    r <- terra::rast(file.path(epath, ENV_FILE))
    ra <- terra::aggregate(r, FACT)
    ra_crp <- terra::crop(ra, RANGE)
    rpts <- terra::extract(ra, RANGE, xy = TRUE)
    if (PTS) {
        return(rpts)
    } else {
        return(ra_crp)
    }
}


#' Extract environmental variables within a grid
#' @description For a given expert range, the extent is calculated and then
#' water is masked out using the world vector. For given grid specs, creates a
#' grid around the extent. Once the grid is created, the environmental
#' variables are extracted.
#' Note: epath is the filepath to the environmental variables, this needs to
#' be in the environment before the function is run
#'
#' @param ENV_FILE list with names of the environmental variables (tiff files)
#' @param FACT factor to aggregate by (default is 10)
#' @param WORLD a world vector (from rnaturaldata)
#' @param RANGE species expert range as a terra vect file
#' @param SEQ length of grid this will make a SEQ x SEQ grid
#' @param RPTS (logical) if TRUE return values else return raster file
#'
#' @return Depending on the value of PTS, either a dataframe or a raster
#'
#' @importFrom terra rast
#' @importFrom terra aggregate
#' @importFrom terra mask
#' @importFrom terra ext
#' @importFrom terra vect
#' @importFrom terra crop
#' @importFrom terra extract
#'
#' @export

extractEnvInExtGrid <- function(ENV_FILE, FACT, WORLD = world, RANGE, SEQ, RPTS) {
    r <- terra::rast(file.path(epath, ENV_FILE))
    ra <- terra::aggregate(r, fact = FACT)
    tmp <- terra::mask(ra, WORLD)
    ex1 <- getExtent(RANGE)
    ex2 <- terra::ext(ex1)
    tmp2 <- terra::crop(tmp, ex2)

    xs <- seq(ex1[1], ex1[2], length.out = SEQ)
    ys <- seq(ex1[3], ex1[4], length.out = SEQ)
    grid <- expand.grid(x = xs, y = ys)
    grid_vect <- terra::vect(grid, geom = c("x", "y"))

    rpts <- terra::extract(tmp2, grid_vect, xy = TRUE)
    rpts <- rpts[complete.cases(rpts), ]

    if (RPTS) {
        return(rpts)
    } else {
        return(tmp2)
    }
}

#' Get extent of a given coordinate of data frame
#' @description This is a little toy function whereas getExtentDF should be
#' used more frequently as it's the "safer" function
#' @param DF dataframe containing coorinates which columns "x" and "y"
#'
#' @return a vector defining extent of coordinates with elements xmin, xmax,
#' ymin, ymax
#'
#' @export

getExtent <- function(DF) {
    xmin <- min(DF$x)
    xmax <- max(DF$x)
    ymin <- min(DF$y)
    ymax <- max(DF$y)
    w <- c(xmin, xmax, ymin, ymax)
    names(w) <- c("xmin", "xmax", "ymin", "ymax")
    return(w)
}

#' Get extent from a data frame
#' @description for a given dataframe containing coordinates, get the extent
#' This is essentially written as a filler function for raster::ext(data.frame)
#' terra's extent function doesn't appear to work for data.frame. This function
#' does that.
#'
#' @param DF data.frame containing coordinates
#' @param NAMES The coordinate columns can be named either c("x", "y") or
#' c("lat", "lon")
#'
#' @return a vector defining extent of coordinates with elements xmin, xmax,
#' ymin, ymax
#'
#' @export

getExtentDf <- function(DF, NAMES) {
    DF2 <- DF[, NAMES]
    # LAT = y, LON = x
    if (all(c("lat", "lon") %in% NAMES)) {
        xmin <- min(DF2[, "lon"])
        xmax <- max(DF2[, "lon"])
        ymin <- min(DF2[, "lat"])
        ymax <- max(DF2[, "lat"])
    }

    if (all(c("x", "y") %in% NAMES)) {
        xmin <- min(DF2[, "x"])
        xmax <- max(DF2[, "x"])
        ymin <- min(DF2[, "y"])
        ymax <- max(DF2[, "y"])
    }
    ext <- c(xmin, xmax, ymin, ymax)
    names(ext) <- c("xmin", "xmax", "ymin", "ymax")

    return(ext)
}

#' Match coordinates between two dataframes
#' @description Takes two dataframes and matches the coordinates in them
#'
#' @param x reference set of coordinates
#' @param match set of coordinates to be matched to x
#'
#' @return An index for match that match x
#'
#' @export

matchCood <- function(x, match) {
    # we want to know which coordinates in x match
    # coordinates in match
    if (class(x)[1] != "data.frame") x <- as.data.frame(x)
    if (class(match)[1] != "data.frame") match <- as.data.frame(match)
    x$id <- paste0(x$x, "_", x$y)
    match$id <- paste0(match$x, "_", match$y)

    ii <- which(match$id %in% x$id)
    return(ii)
} # Returns ratio for the coordinate system
gimmeRatios <- function(YFSP, YCOOD, REF_RAST) {
    # Parition the YFSP
    # First quick check
    if (length(YFSP) != nrow(YCOOD)) {
        error("Something is wrong - Y and COOD are not equal dimensions")
    }
    ind1 <- which(YFSP == 1)
    ind0 <- which(YFSP == 0)

    # Dealing with absences
    # Create a raster with counts of absences coordinates/cell
    c0 <- YCOOD[ind0, ]
    c0_vect <- terra::vect(c0[, c(2, 1)], geom = c("lon", "lat"))
    er0 <- terra::rasterize(c0_vect, REF_RAST, fun = sum)
    v0 <- terra::values(er0)
    inv0 <- which(is.na(v0))
    v0[inv0] <- 0
    terra::values(er0) <- v0

    # Dealing with presences
    # Create a raster with counts of presences coordinates/cell
    c1 <- YCOOD[ind1, ]
    c1_vect <- terra::vect(c1[, c(2, 1)], geom = c("lon", "lat"))
    er1 <- terra::rasterize(c1_vect, REF_RAST, fun = sum)
    v1 <- terra::values(er1)
    inv1 <- which(is.na(v1))
    v1[inv1] <- 0
    terra::values(er1) <- v1

    # Calculate ratio
    err <- (er1 / er0)
    vr <- terra::values(err)
    invr <- which(is.infinite(vr))
    vr[invr] <- 1
    invr0 <- which(vr == 0)
    vr[invr0] <- NA
    terra::values(err) <- vr

    # dfr <- terra::as.data.frame(err, xy = TRUE)
    p <- terra::as.points(err)
    data <- cbind(terra::crds(p), terra::values(p))

    return(data)
}




#' Make a grid for a given set of coordinates
#' @description Take a dataframe containing coordinates and convert to a
#' SpatialPoints class (sp packages -- need to update possibly) and then
#' use that to make a grid for a specified cellsize
#'
#' @param PRJ projection string, default WGS84
#' @param COORDS coordinate dataframe
#' @param CELLSIZE specify grid cell size
#'
#' @return a terra vector grid
#'
#' @importFrom sp SpatialPoints
#' @importFrom sp CRS
#' @importFrom sf st_as_sf
#' @importFrom sf st_crs
#' @importFrom sf st_bbox
#' @importFrom sf st_as_sfc
#' @importFrom sf st_make_grid
#' @importFrom terra vect
#'
#' @export

makeGrid <- function(PRJ = "+proj=longlat +datum=WGS84", COORDS,
                     CELLSIZE = .2) {
    # A little housekeeping - getting coordinate dataframe ready
    coords <- sp::SpatialPoints(COORDS, proj4string = sp::CRS(PRJ))
    data_sf <- sf::st_as_sf(coords,
        coords = c("long", "lat"),
        crs = sf::st_crs(PRJ)
    )

    # Create a grid
    grid <- data_sf |>
        sf::st_bbox() |>
        sf::st_as_sfc() |>
        sf::st_make_grid(
            cellsize = c(CELLSIZE, CELLSIZE),
            crs = PRJ,
            square = TRUE
        ) |>
        sf::st_as_sf()
    grid <- terra::vect(grid)

    return(grid)
}

#' Distance calculation between two sets of coordinates
#' @description pairwise distance calculation for specifically two sets of
#' coordinates and return the minimum distance between them for example, one
#' set of centroids from the region and one set of centroids for the range
#'
#' @param SFP absolute Script File Path
#' @param CEN_GRID centroids (or points) from the grid/region
#' @param CEN_RANGE centroids (or points)
#'
#' @return Return a matrix of minimum distance of length nrow(CEN_GRID)
#'
#' @importFrom Rcpp sourceCpp
#'
#' @export

distanceCalc <- function(SFP = file.path(root, "scripts/"),
                         CEN_GRID, CEN_RANGE) {
    # Source the cpp function
    Rcpp::sourceCpp(paste0(SFP, "haversine_vec.cpp"))
    # Set up for loop
    ncen <- nrow(CEN_GRID)
    dist <- matrix(0, nrow = ncen)
    ncen_range <- nrow(CEN_RANGE)
    # Loop
    for (i in seq_len(ncen)) {
        region_centroids <- data.frame(
            x = rep(CEN_GRID[i, "x"], ncen_range),
            y = rep(CEN_GRID[i, "y"], ncen_range)
        )
        # Now calculate distances
        tmp <- haversineVec(
            lat1 = region_centroids[, "y"], lon1 = region_centroids[, "x"],
            lat2 = CEN_RANGE[, "y"], lon2 = CEN_RANGE[, "x"]
        )
        # Get minimum distance
        dist[i] <- min(tmp)
    }

    return(dist)
}

#' Logistic function
#' @description Pretty self explanatory
#' X1 = L / (1 + exp(-K x DIST - X0))
#'
#' @param X0 Mid distance
#' @param L Upper limit
#' @param K Decay parameter
#' @param DIST Distances
#'
#' @return Output of the logisitic function same length as DIST
#'
#' @export

logisticFunction <- function(X0, L, K, DIST) {
    x <- (L / (1 + exp(-K * (DIST - X0))))
    return(x)
}

# pred_functions ------
# Functions for prediction

# Maximize sensitivity
#
# For a given number of thresholds, maximize sensitivity values
#
# @param XCONT Predictions (probability of presence i.e. continuous on [0,1])
# @param ACTUAL Observed
# @param ALL Return entire dataframe (TRUE) or just the maximized value
#
# @return either a dataframe or the maximized sensitivity value
#
# @export

# maxSens <- function(XCONT, ACTUAL, ALL = T) {
#
#   # Set up some thresholds
#   thresholds <- seq(0.1, 0.9, by = 0.1)
#   df <- data.frame(thresholds = thresholds,
#                    sens = rep(0, length(thresholds)))
#
#   # Test sensitivity per threshold
#
#   for(i in 1:length(thresholds)) {
#     xbin <- ifelse(XCONT > quantile(XCONT, thresholds[i]), 1, 0)
#     if(sum(xbin) != 0) {
#       conf_matrix <- table(ACTUAL, xbin)
#       df[i, "sens"] <- caret::sensitivity(conf_matrix)
#     }
#   }
#
#   if(!ALL) {
#     r <- which(df[, 'sens'] == max(df[, "sens"]))
#     return(df[r,])
#   } else {
#     return(df)
#   }
#
# }


# Maximize specificity
#
# For a given number of thresholds, maximize specificity values
#
# @param XCONT Predictions (probability of presence i.e. continuous on [0,1])
# @param ACTUAL Observed
# @param ALL Return entire dataframe (TRUE) or just the maximized value
#
# @return either a dataframe or the maximized specificity value
#
# @export


# maxSpec <- function(XCONT, ACTUAL, ALL = T) {
#
#   # Set up some thresholds
#   thresholds <- seq(0.1, 0.9, by = 0.1)
#   df <- data.frame(thresholds = thresholds,
#                    spec = rep(0, length(thresholds)))
#
#   # Test sensitivity per threshold
#
#   for(i in 1:length(thresholds)) {
#     xbin <- ifelse(XONT > quantile(XCONT, thresholds[i]), 1, 0)
#     if(sum(xbin) != 0) {
#       conf_matrix <- table(ACTUAL, xbin)
#       df[i, "spec"] <- caret::specificity(conf_matrix)
#     }
#   }
#
#   if(!ALL) {
#     r <- which(df[, 'spec'] == max(df[, "spec"]))
#     return(df[r,])
#   } else {
#     return(df)
#   }
#
# }

#' Make a Beta Matrix
#'
#' This matrix holds all the posterior means from the beta posteriors and has
#' dimensions K x J
#'
#' @param BETA list of beta means with length K
#'
#' @return beta matrix of dimensions K X J
#'
#' @export

makeBetaMat <- function(BETA) {
    k <- length(BETA)
    x <- lapply(BETA, colMeans)
    j <- length(x[[1]])

    bmat <- do.call(rbind, x)

    return(bmat)
}

#' Get conditionally predicted files
#'
#' Read in and organize species-specific files into a matrix
#'
#' @param GENUS which is the focal genus
#' @param VAR options are ENV (for phylogenetically correlated vars),
#' EFFORT (for non-phylogenetically correlated vars)
#' @param EXP_FLAG which experiment are you running e.g. ALL, NOCC, NOPHYLO,
#' DEF5, etc
#' @param DIR the result directory
#' @param SPNO Number of species
#' @param NEFF Number effort variables
#' @param NENV Number of environmental variables
#' @param SP Species names
#' @param REP rep number
#' @param DEFLEV the deficiency level [nodef, y100, y25 etc]. Except for nodef
#' the deflevel should be followed by y and the function assumes in the file
#' name, the deficiency level is surrounded by punctuation e.g. _y50_
#'
#' @return conditional predicted beta matrix of either effort or env variables
#' for specified experiment and genus
#'
#' @export

getCondPred <- function(GENUS = genus, VAR, EXP_FLAG,
                        DIR = file.path(res.directory, "condPred_res"),
                        NEFF = E, NENV = Konly, SPNO = J, SP = sp,
                        REP = rep, DEFLEV = deflev) {
    # Get files and clean - remove sd files
    path <- file.path(DIR, GENUS)
    files <- dir(path)

    clean <- grep("sd", files)
    if (length(clean) > 0) files <- files[-clean]

    fin <- grep(VAR, files)
    files <- files[fin]

    fin2 <- grep(EXP_FLAG, files)
    files <- files[fin2]

    fin3 <- grep(paste0("rep", REP), files)
    files <- files[fin3]

    rexp <- paste0("[[:punct:]]", DEFLEV, "[[:punct:]]", collapse = "")
    fin4 <- grep(rexp, files)
    files <- files[fin4]

    # Check - this should be true
    if (length(files) != SPNO) {
        warning("Something is wrong, length of files != number of species")
    }

    if (grepl("ENV", VAR)) R <- NENV
    if (grepl("EFF", VAR)) R <- NEFF
    # if(!(VAR %in% c("ENV", "EFF"))) stop("Don't recognize VAR input")

    rin <- data.frame("no" = seq_along(files), "files" = files)
    store <- apply(rin, 1, function(x) {
        w <- read.csv(file.path(path, x["files"]))
        return(w)
    })
    x <- lapply(store, colMeans)
    xmat <- matrix(unlist(x), ncol = SPNO, nrow = R)
    colnames(xmat) <- SP

    return(xmat)
}

#' (Pred)iction at given (loc)ation function
#'
#' For a given matrix of Xs, betas and variances, function returns
#' predicted presence probability in the given locations or predicted
#' presence/absence. It also calculates the uncertainty at every location
#'
#' @param X the covariate matrix of dimension N x K
#' @param BMAT the beta matrix with posterior beta means of dimension K x J
#' @param YDEVS y standard deviations posterior means
#' @param BIN return predicted presence/absence (TRUE) or predicted
#' probability of presence (FALSE)
#' @param SP species names
#'
#' @return a list of either predicted presence/absence or probability presence
#' and standard deviation of predictions (at every location)
#'
#' @importFrom MASS mvrnorm
#' @importFrom mvtnorm pmvnorm
#'
#' @export

predLoc <- function(X, BMAT, YDEVS, BIN, SP = sp) {
    xbeta <- X %*% BMAT

    j <- length(YDEVS)
    cov <- matrix(0, j, j)
    diag(cov) <- YDEVS

    z <- list()
    zsd <- list()
    un <- list()

    for (i in seq_len(nrow(xbeta))) {
        if (BIN) {
            # If binary, return predicted presence/absence
            z[[i]] <- MASS::mvrnorm(n = 100, mu = xbeta[i, ], Sigma = cov)
            z[[i]] <- colMeans(z[[i]])
            z[[i]] <- ifelse(z[[i]] > 0, 1, 0)
        } else {
            # Else return predicted probability of presence
            x <- MASS::mvrnorm(n = 1, mu = xbeta[i, ], Sigma = cov)
            tmp <- c()

            for (ii in seq_along(x)) {
                tmp[ii] <- mvtnorm::pmvnorm(x[ii], mean = xbeta[i, ], sigma = cov)[1]
            }
            z[[i]] <- tmp
        }

        # Get predicted uncertainty
        un[[i]] <- MASS::mvrnorm(n = 100, mu = xbeta[i, ], Sigma = cov)
        zsd[[i]] <- apply(un[[i]], 2, sd)
    }

    zmat <- do.call(rbind, z)
    zsd <- do.call(rbind, zsd)

    colnames(zmat) <- colnames(zsd) <- SP
    return(list("zmat" = zmat, "zsd" = zsd))
}

#' Compute AUC values for given predicted and observed values
#'
#' @param Y matrix of observed presence/absence values (N x J)
#' @param PRED matrix of predicted presence/absence values (N x J)
#' @param PRINT write out print statements - this is kind of annoying
#'
#' @return return a vector of AUC values for given number of species
#'
#' @importFrom pROC auc
#'
#' @export

computeAUC <- function(Y = y, PRED, PRINT = FALSE) {
    nsp <- ncol(Y)
    x <- matrix(0, nrow = 1, ncol = nsp)

    # Check point
    if (all(colnames(PRED) != colnames(Y))) {
        warning("Columns of prediction might not match up with columns of y")
    }

    for (j in seq_len(nsp)) {
        x[1, j] <- suppressMessages(pROC::auc(Y[, j], PRED[, j]))
    }


    colnames(x) <- colnames(Y)
    return(x)
}


#' Compute sensitivity (hit rate)
#'
#' Computes the proportion of correctly predicted presences
#'
#' @param Y observed data
#' @param PRED predicted data
#'
#' @return sensitivity value for all ys
#'
#' @importFrom caret sensitivity
#'
#' @export

computeSens <- function(Y = y, PRED) {
    nsp <- ncol(Y)
    x <- matrix(0, nrow = 1, ncol = nsp)

    for (j in 1:nsp) {
        conf_matrix <- table(PRED[, j], Y[, j])
        x[, j] <- caret::sensitivity(conf_matrix)
    }

    return(x)
}


#' Compute specificity
#'
#' Computes the proportion of correctly predicted absences
#'
#' @param Y observed data
#' @param PRED predicted data
#'
#' @return specificity values for all ys
#'
#' @importFrom caret specificity
#'
#' @export

computeSpecs <- function(Y = y, PRED) {
    nsp <- ncol(Y)
    x <- matrix(0, nrow = 1, ncol = nsp)

    for (j in 1:nsp) {
        conf_matrix <- table(PRED[, j], Y[, j])
        x[, j] <- caret::specificity(conf_matrix)
    }

    return(x)
}


#' Compute rate of false positives
#'
#' Computes the proportion of false positives
#'
#' @param Y2 observed data
#' @param PRED2 predicted data
#'
#' @return false positive proportions values
#'
#' @importFrom caret specificity
#'
#' @export

computeFalsePos <- function(Y2, PRED2) {
    specs <- computeSpecs(Y2, PRED2)
    x <- 1 - specs

    return(x)
}


#' Make a labelled confusion matrix

#' Creates a confusion matrix - can't be fed directly into caret
#' functions because of labels
#'
#' @param Y observed data
#' @param PRED predicted data
#'
#' @return list of labelled confusion matrix
#'
#' @export

makeConfusionMatrix <- function(Y = y, PRED) {
    xx <- list()
    nsp <- ncol(Y)

    for (i in 1:nsp) {
        x <- table(PRED[, i], Y[, i])
        colnames(x) <- c("actual0", "actual1")
        rownames(x) <- c("pred0", "pred1")
        xx[[i]] <- x
    }


    return(xx)
}


#' Compute frequency bias
#'
#' Compares elements of confusion matrix to assess model under/over prediction
#'
#' If A = true presence, B = false presence, C = false absence, D = true absence
#' model formula is (A + B)/(A + C)
#'
#' 0.5 is neither over nor under. 0.5 + is over, 0.5 - is under prediction
#'
#' @param Y observed data
#' @param PRED predicted data
#'
#' @return bias value
#'
#' @export

computeFreqBias <- function(CONF_MAT) {
    bias <- (CONF_MAT[4] + CONF_MAT[2]) / (CONF_MAT[4] + CONF_MAT[3])

    return(bias)
}

#' Compute prevalence
#'
#' Measures the ratio of presences (both predicted and omitted) to all cells
#' and hence expresses how common within the study area a species is
#'
#' If A = true presence, B = false presence, C = false absence, D = true absence
#' model formula is (A + C)/(A + B + C + D)
#'
#' 1 is prevalent, 0 is not
#'
#' @param CONF_MAT confusion matrix
#'
#' @return prevalence value
#'
#' @export

computePrevalence <- function(CONF_MAT) {
    p <- (CONF_MAT[4] + CONF_MAT[3]) / (sum(CONF_MAT))
    return(p)
}


#' Compute True Skill Score (TSS)
#'
#' Measured as Sensitivity - (1 - Specificity)
#'
#' @param Y2 observed y
#' @param PRED2 predicted y
#'
#' @return tss value
#'
#' @export


computeTSS <- function(Y2, PRED2) {
    sens <- computeSens(Y2, PRED2)
    fp <- computeFalsePos(Y2, PRED2)

    tss <- sens - fp

    return(tss)
}

#' Compute Odd Ratio Skill Score (ORSS)
#'
#' If A = true presence, B = false presence, C = false absence, D = true absence
#' model formula is (AD - BC)/(AD + BC)
#'
#' @param CONF_MAT confusion matrix
#'
#' @return orss value
#'
#' @export


computeORSS <- function(CONF_MAT) {
    t1 <- CONF_MAT[4] * CONF_MAT[1]
    t2 <- CONF_MAT[2] * CONF_MAT[3]

    orss <- (t1 - t2) / (t1 + t2)

    return(orss)
}

#' Store a prediction
#'
#' Make and store a prediction (file structure assumed and used within
#' the function)
#'
#' @param x1 X covariate matrix
#' @param bmat3 A beta matrix
#' @param ydevs1 a vector length J of observation error
#' @param binary (logical) TRUE if binary, FALSE if continuous
#' @param fsp1 focal species name
#' @param sp1 character vector length J with names of species
#'
#' @return Nothing, just
#' @export

storePred <- function(x1 = x, bmat3 = bmat2, ydevs1 = ydevs,
                      binary = TRUE, fsp1 = fsp, sp1 = sp, wpath, rep1 = rep) {
    # Make a prediction
    predA <- predLoc(x1, bmat3, ydevs1, TRUE, sp1)
    pbin <- predA$zmat
    pzd <- predA$zsd

    # Store a prediction
    ppath <- file.path(wpath, "spatial_pred")
    write.csv(pbin, file = file.path(ppath, paste0(
        fsp1, "_",
        exp_flag, "rep", rep1, "_presence_absence.csv"
    )))
    write.csv(pbin, file = file.path(ppath, paste0(
        fsp1, "_",
        exp_flag, "rep", rep1, "_uncertainty.csv"
    )))
}

#' Evaluate a prediction and store AUC and other performance metrics
#' This function is basically a wrapper that takes in x and y inputs,
#' uses the pred function to create a prediction and then uses auc (and other)
#' metric functions to produce evaluation metrics.
#'
#' @param x1 the x (covariate) data
#' @param y1 the y (observed) data
#' @param bmat3 is the beta matrix with dimensions K x J
#' @param ydevs1 is a vector of observation error of length J
#' @param bin (logical) TRUE for binary FALSE for continuous data
#' @param sp1 is a vector of all the species names
#' @param fsp1 is a string with the name of the focal species
#' @param wpath is the write out path
#' @param insamp (logical) TRUE for insample prediction, FALSE for out of sample
#' prediction
#' @param writeout (logical) TRUE for write out FALSE o/w
#' @param iter numeric value for number of iterations for AUC and other
#' metric calculation
#' @param rep1 the rep number
#' @param deflev1 the deficiency level [nodef, y100, y25 etc]. Except for nodef
#' the deflevel should be followed by y and the function assumes in the file
#' name, the deficiency level is surrounded by punctuation e.g. _y50_
#'
#' @return if writeout, the files are written out o/w nothing is returned
#'
#' @export

evaluatePred <- function(x1, y1, bmat3 = bmat2, ydevs1 = ydevs,
                         bin = TRUE, sp1 = sp, fsp1 = fsp,
                         wpath = file.path(path1), insamp,
                         writeout = TRUE, iter = 100, rep1 = rep,
                         deflev1 = deflev) {
    J <- ncol(bmat3)
    # Loop of aucs
    iter <- 100
    auc <- matrix(0, ncol = J, nrow = iter)
    sens <- matrix(0, ncol = J, nrow = iter)
    bias <- matrix(0, ncol = J, nrow = iter)
    tss <- matrix(0, ncol = J, nrow = iter)

    for (i in 1:iter) {
        # prediction
        predA <- predLoc(x1, bmat3, ydevs1, bin, sp1)
        pbin <- predA$zmat

        # Won't work without only 0s
        if (any(colSums(pbin) == 0)) {
            print("Here -- 0")
            ind <- which(colSums(pbin) == 0)
            pbin[sample(seq_len(nrow(pbin)), size = 1), ind] <- 1
        }

        # Won't work without only 1s
        if (any(colSums(pbin) == nrow(pbin))) {
            print("Here --1")
            ind <- which(colSums(pbin) == nrow(pbin))
            pbin[sample(seq_len(nrow(pbin)), 1), ind] <- 0
        }

        # print(i)
        # AUCs and other metrics
        auc[i, ] <- computeAUC(Y = y1, PRED = pbin, PRINT = FALSE)
        # sens[i, ] <- computeSens(Y = y1, PRED = pbin)
        # tss[i, ] <- computeTSS(Y2 = y1, PRED2 = pbin)
        # conf_matrix <- makeConfusionMatrix(Y = y, PRED = pbin)
        # bias[i, ] <- c(unlist(lapply(conf_matrix, computeFreqBias)))
    }

    print(colMeans(auc))
    # colMeans(sens)
    # colMeans(tss)
    # colMeans(bias)
    if (insamp) {
        trail <- "INSAMPLE"
    } else {
        trail <- "OUTSAMPLE"
    }

    if (writeout) {
        write.csv(auc,
            file = file.path(wpath, "auc", paste0(
                fsp1, "_", exp_flag,
                "_", deflev, "_", trail, "_", "rep", rep1, "_auc.csv"
            ))
        )
        # write.csv(sens, file.path(wpath, "sens",
        #                          paste0(focal_sp, "_", exp_flag, "_sens.csv")))
        # write.csv(tss, file.path(wpath, "tss",
        #                          paste0(focal_sp, "_", exp_flag, "_tss.csv")))
        # write.csv(bias, file.path(wpath, "bias",
        #                          paste0(focal_sp, "_", exp_flag, "_bias.csv")))
    }
}


# raw_data_functions -----

# Raw data wrangling functions


#' Organize data as cood-y
#' @description Clean up the raw data and return the
#' columns required and names formatted
#'
#' @param(DPATH data path
#' @param name name of the dataset
#' @param columns name of columns apart from species columns to be returned
#'
#' @return dataset with selected colums
#'
#' @export

organizeY <- function(DPATH, name, columns = c(
                          "lat", "lon", "observation_date",
                          "num.observers", "duration", "distance"
                      )) {
    if (name == "hummingbird_sa") {
        print("1. Loading data...")
        raw <- read.csv(file.path(DPATH, paste0(name, ".csv")))

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

        return(raw[, c(columns, spnm)]) # nolint
    }
}

#' Check that phylogenetic tree names and occurrence species names are the
#' same
#'
#' @param tree is the phylogenetic tree, object of class phylo
#' @param y matrix or dataframe with occurrences, organized as site x species
#'
#' @return a list with elements y and tree where they are the same classes as
#' the input but now with the same species between them
#'
#' @export
#'

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
    if (all(ysp %in% tsp) && all(tsp %in% ysp)) {
        print("All species names match")
    } else {
        print("Something went wrong - debug with y <- Y; tree <- TREE")
    }

    return(list(y = y, tree = tree))
}

#' Organize data as cood-y
#'
#' @param(DPATH absolute path to where the data resides
#' @param NAME name of dataset to be read in
#'
#' @return a dataframe with lat, lon, obs date, effort and species columns
#' formatted by replacing "." with "_"
#'
#' @export
#'

organizeY <- function(DPATH, NAME) {
    if (NAME == "hummingbird_sa") {
        print("1. Loading data...")
        raw <- read.csv(file.path(DPATH, paste0(NAME, ".csv")))

        print("2. Organizing data...")
        # Replace "." in sp names with "_"
        # Species start at column 8
        # TODO: this sucks and will definitely f*** things up later
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

#' Check phylo sp and occ sp are the same
#'
#' @param TREE phylogenetic tree (object from class phylo)
#' @param Y dataset with all the species records
#'
#' @return list with elements "tree" and "y" that have harmonized names
#'
#' @importFrom ape drop.tip
#'
#' @export

equalizeSP <- function(TREE, Y) {
    # y needs to be just species records
    tsp <- TREE$tip.label
    ysp <- colnames(Y)

    # Which y sp NOT IN t sp
    ky <- which(!(ysp %in% tsp))
    if (!identical(ky, integer(0))) Y <- Y[, -c(ky)]

    # Which t sp NOT IN y sp
    ysp <- colnames(Y)
    kt <- which(!(tsp %in% ysp))
    if (!identical(kt, integer(0))) TREE <- ape::drop.tip(TREE, kt)

    # Final check
    tsp <- TREE$tip.label
    if (all(ysp %in% tsp) && all(tsp %in% ysp)) {
        print("All species names match")
    } else {
        print("Something went wrong - debug with y <- Y; tree <- TREE")
    }

    return(list(y = Y, tree = TREE))
}


#' Organize X covariates
#'
#' @param EPATH path to environmental rasters
#' @param ENV_FILES vector with string of all the names of env rasters
#' @param ENV_CRS crs string
#' @param COOD cood dataframes with column names ("lon", "lat")
#'
#' @export
#'
#' @importFrom raster extent
#' @importFrom raster stack
#' @importFrom raster compareCRS
#' @importFrom raster raster
#' @importFrom raster crop
#' @importFrom raster projectRaster
#' @importFrom raster crs
#' @importFrom raster addLayer
#' @importFrom raster extract

# TODO: replace raster with terra!

# organizeCovar <- function(EPATH, ENV_FILES, ENV_CRS, COOD) {
#     env.paths <- file.path(EPATH, ENV_FILES)
#     colnames(COOD) <- c("y", "x")
#     COOD <- COOD[, c("x", "y")]
#     extent <- terra::ext(COOD)

#     env <- raster::stack()
#     for (env.path in env.paths) {
#         print(paste0("Loading ", env.path, " ..."))
#         env.layer <- raster::raster(env.path)

#         # If needed, reproject the env layer
#         if (!raster::compareCRS(raster::crs(env.layer), ENV_CRS)) {
#             print(paste0(" Reprojecting..."))
#             env.layer <- raster::projectRaster(env.layer, ENV_CRS)
#         }

#         if (!is.null(extent)) {
#             env.crop <- raster::crop(env.layer, extent)
#         }
#         env <- raster::addLayer(env, env.crop)
#     }

#     cood.mat <- do.call(cbind, COOD)

#     annotations <- raster::extract(env, cood.mat,
#         method = "bilinear"
#     )
#     return(annotations)
# }

#' Function for creating a reference coord set
#'
#' @param REF_RAST reference raster for a blank dataframe with coordinates
#'
#' @importFrom rnaturalearth ne_countries
#' @importFrom terra vect
#' @importFrom terra mask
#' @importFrom terra as.data.frame
#'
#' @export

makeRefCoods <- function(REF_RAST, AREA = area) {
    REF_RAST[] <- 1
    area_vect <- terra::vect(AREA)
    # Mask
    REF_RAST <- terra::mask(REF_RAST, area_vect)
    w <- terra::as.data.frame(REF_RAST, xy = TRUE)
    # Return reference set of coordinates
    return(w)
}

#' Calculate the ratios of presences to absences in every grid cell
#'
#' @param YFSP data for y focal species
#' @param YCOOD data for y focal species coordinates
#' @param REF_RAST reference raster
#'
#' @importFrom terra cellFromXY
#' @importFrom terra as.data.frame
#'
#' @export

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

# Sampling functions ------

# Sampling functions

#' Get extent and widen around the presences of the focal species
#'
#' @param Y2P coordinates and y
#'
#' @return an object of class extent indicating modelling extent
#'
#' @importFrom raster extent
#'
#' @export


widenExtent <- function(Y2C) {
    y2 <- Y2C[, -c(1, 2)]
    # Create presence only dataset
    i <- which(rowSums(y2) > 0)
    y2p <- Y2C[i, ]

    # Crop the whole dataset to a smaller region
    ex <- terra::ext(y2p)
    widex <- ex + 1

    return(widex)
}

#' Sub sample the full EBird hummingbird dataset to get a better
#' proportion of presences and absences
#'
#' @param COOD coordinates from the full dataset
#' @param Y2 focal species/clade
#' @param X environmental data annotating full dataset
#' @param FSP (Optional) if not null, i.e. if a species index is provided,
#' the dataset is subsampled trying to optimize the dataset for that species
#' @param THIN_ALL TRUE/FALSE - thin entire clade to customize for the focal
#' species ?
#' (yes/TRUE)
#' @param PERC presence-absence ratio (50% is the default)
#'
#' @return a list of subsampled, x, y and coods
#'
#' @importFrom raster extent
#'
#' @export


subSample <- function(COOD, Y2, X, FSPIND = NULL,
                      THIN_ALL = FALSE, PERC = 50) {
    y2c <- cbind(COOD, Y2)
    colnames(y2c)[c(1, 2)] <- c("y", "x")
    c2 <- y2c[, c("x", "y")]
    i <- which(rowSums(Y2) > 0)
    y2p <- y2c[i, ]
    widex <- widenExtent(y2p)
    cood_crop <- which(c2$x > widex[1] & c2$x < widex[2] & c2$y >
        widex[3] & c2$y < widex[4])
    ycrp <- y2c[cood_crop, ]
    np <- colSums(ycrp[, 3:ncol(ycrp)])
    nt <- np * (100 / PERC)
    na <- nt - np
    pool <- ycrp[, c(3:ncol(ycrp))]
    a_cood <- ycrp[which(rowSums(pool) == 0), ]
    if (is.null(FSPIND)) {
        na_size <- mean(na)
    } else {
        na_size <- na[FSPIND]
    }
    ind <- sample(seq_along(a_cood), na_size)
    a_sel <- a_cood[ind, ]
    if (THIN_ALL) {
        if (is.null(FSPIND)) {
            warning("Focal species index can't be null")
        }
        xy <- y2p[, c("x", "y")]
        y2 <- y2p[, -c(1, 2)]
        pind <- which(y2[, FSPIND] == 1)
        aind <- which(y2[, FSPIND] == 0)
        op <- data.frame(
            ind = which(rowSums(y2[aind, ]) > 0),
            sums = rowSums(y2[aind, ])
        )
        op <- op[order(op$sums, decreasing = TRUE), ]
        if (nrow(op) > na_size) {
            ww <- op$ind[1:na_size]
            yfin <- y2p[c(pind, ww), ]
        } else {
            m <- na_size - nrow(op)
            ind <- sample(seq_along(a_cood), m)
            a_sel <- a_cood[ind, ]
            yfin <- rbind(y2p, a_sel)
        }
    } else {
        yfin <- rbind(y2p, a_sel)
    }
    cood3 <- yfin[, c("y", "x")]
    cs <- paste0(COOD$y, "_", COOD$x)
    cs3 <- paste0(cood3$y, "_", cood3$x)
    ind <- match(cs3, cs)
    x2 <- X[ind, ]
    y3 <- yfin[, c(3:ncol(yfin))]
    return(list(x = x2, y = y3, cood = cood3))
}


bgSample <- function(COOD, Y2, X, FSPIND = NULL, THIN_ALL = FALSE,
                     PERC = 50) {
    y2c <- cbind(COOD, Y2)
    colnames(y2c)[c(1, 2)] <- c("y", "x")
    c2 <- y2c[, c("x", "y")]

    # Create presence only dataset
    i <- which(rowSums(Y2) > 0)
    y2p <- y2c[i, ]

    # Crop the whole dataset to a smaller region

    widex <- makeDomain(y2p)
    cood_crop <- which(c2$x > widex[1] & c2$x < widex[2] &
        c2$y > widex[3] & c2$y < widex[4])

    ycrp <- y2c[cood_crop, ]

    # Calculate proportions of absences needed
    np <- colSums(ycrp[, 3:ncol(ycrp)])

    # Back calculate how many absences we need for given p/a ratio
    nt <- np * (100 / PERC)
    na <- nt - np

    pool <- ycrp[, c(3:ncol(ycrp))]
}



#' Create a data-deficient dataset for a given focal species
#'
#' @param FSPIND focal species index
#' @param FSP focal species name
#' @param Y occurrence dataset without coods
#' @param COOD coordinates
#' @param X environmental variables
#' @param DEF_LEV level of deficiency - "F" for no deficiency
#'
#' @return list of data-deficient y, x, and coods
#'
#' @export

create_deficiency <- function(FSPIND = fspind, FSP = fsp,
                              Y = y3, X = x2, COOD = cood2,
                              DEF_LEV) {
    npres <- sum(Y[, FSPIND])

    pres_ind <- which(Y[, FSPIND] == 1)
    size <- length(pres_ind) - DEF_LEV
    y5_ind <- sample(pres_ind, size)

    y5 <- Y[-y5_ind, ]
    x5 <- X[-y5_ind, ]
    cood5 <- COOD[-y5_ind, ]

    return(list(
        "x" = x5,
        "y" = y5,
        "cood" = cood5
    ))
}


#' Create a domain as union of buffered presence points and expert range map
#'
#' @param Y2 matrix or dataframe containing the occurrence points (must include
#' the focal species data)
#' @param FSPIND focal species index (specifically which column in the Y2 data
#' is the focal species)
#' @param COOD matrix or dataframe containing coordinates of sampled points -
#' must be the same number of rows as Y@ and have column names "y" and "x"
#' for lat-lon
#'
#' @param EXPERT_PATH path to expert ranges
#' @param WRITE_DOMAIN logical TRUE if domains should be written out to domain
#' sub-directory
#' @param BUFFER_EX buffer extent (in meters)
#' @param UTMZONE utm zone (17 for SA)
#' @param PLOT_DOMAIN logical TRUE if want to visualize
#'
#' @import sf
#' @import sp
#'
#' @export
#'


createDomain <- function(Y2, FSPIND, COOD,
                         EXPERT_PATH =
                             "~/phylo-sdms/phyloproj/raw_data/expert_ranges",
                         WRITE_DOMAIN = TRUE,
                         DOMAIN_PATH =
                             "~/phylo-sdms/phyloproj/raw_data/domains",
                         BUFFER_EX = 5e05,
                         UTMZONE = 17,
                         PLOT_DOMAIN = TRUE) {
    # Set domain as the union of 5000m polygons around presence points
    # and expert range map
    library(sf)
    library(sp)
    library(raster)


    # Only need presence points for finding the modelling domain
    k <- which(Y2[, FSPIND] == 1)
    fsp <- colnames(Y2)[FSPIND]

    y3 <- y2[k, FSPIND]
    cood2 <- COOD[k, ]
    yc2 <- data.frame(cbind(cood2, y3))
    colnames(yc2) <- c("y", "x", colnames(y2)[FSPIND])
    pnts_sf <- sf::st_as_sf(yc2, coords = c("x", "y"), crs = sf::st_crs(4326))


    # Read in expert range map shapefile
    r <- sf::read_sf(file.path(EXPERT_PATH, paste0(fsp, ".shp")))
    r <- sf::st_as_sf(r)

    # Code for intersections -

    # UTM Zone 17N roughly
    utmStr <- "+proj=utm +zone=%d +datum=NAD83 +units=m +no_defs +ellps=GRS80"
    crs <- sp::CRS(sprintf(utmStr, UTMZONE))
    pnts_trans <- sf::st_transform(pnts_sf, 4236)
    r_trans <- sf::st_transform(r, 4236)
    # pnts_trans <- pnts_sf %>% mutate(
    # intersection = as.integer(st_intersects( pnts_trans, r_trans)))
    # pnts_trans

    # Code for buffer
    buffers <- sf::st_buffer(pnts_trans, BUFFER_EX) # 5km around each p point

    # Code to take the union with expert ranges
    buffers_trans <- sf::st_transform(buffers, 4236)
    domain <- st_union(r_trans, buffers_trans)

    p <- ggplot() +
        geom_sf(data = domain)
    if (PLOT_DOMAIN) p

    if (WRITE_DOMAIN) {
        save(domain, file = file.path(
            DOMAIN_PATH,
            paste0("domain_", fsp, ".Rdata")
        ))
    }

    return(domain)
}

#' Get Presence IDs
#' @description returns presence ids randomly sampled from a pool given a size
#'
#' @param YC dataframe with "lat", "lon", species observations, and "ids"
#' @param SP species index or name
#' @param SIZE number of presence ids to be sampled and returned - must be
#' less than number of rows in YC
#'
#' @return vector of ids
#'
#' @export

getPresenceIDs <- function(YC, SP, SIZE) {
    x <- YC[, c("lon", "lat", SP, "id")]
    xx <- x[which(x[, SP] == 1), ]
    if (nrow(xx) > SIZE) {
        size1 <- SIZE
    } else {
        size1 <- nrow(xx)
    }
    k <- sample(seq_len(nrow(xx)), size1, replace = FALSE)
    ids <- xx[k, "id"]
    return(ids)
}


#' gimmeRatios
#' @description returns the coordinates and ratio of presence - absences at that
#' location
#'
#' @param YFSP vector with presences/absences (1/0s) for the focal species
#' @param YCOOD dataframe or matrix with coordinates (specifically "lon", "lat")
#' that correspond to the observations in YFSP.
#' @param REF_RAST reference raster that will be used to rasterize the coods
#' and the returned cood + ratio will be from this raster
#'
#' @return data frame with ratio values and coordinates
#'
#' @importFrom terra vect
#' @importFrom terra rasterize
#' @importFrom terra values
#' @importFrom terra as.points
#' @importFrom terra crds
#'
#' @export


gimmeRatios <- function(YFSP, YCOOD, REF_RAST) {
    # Parition the YFSP
    # First quick check
    if (length(YFSP) != nrow(YCOOD)) {
        error("Something is wrong - Y and COOD are not equal dimensions")
    }
    ind1 <- which(YFSP == 1)
    ind0 <- which(YFSP == 0)

    # Dealing with absences
    # Create a raster with counts of absences coordinates/cell
    c0 <- YCOOD[ind0, ]
    c0_vect <- terra::vect(c0[, c(2, 1)], geom = c("lon", "lat"))
    er0 <- terra::rasterize(c0_vect, REF_RAST, fun = sum)
    er0[is.na(er0)] <- 0

    # Dealing with presences
    # Create a raster with counts of presences coordinates/cell
    c1 <- YCOOD[ind1, ]
    c1_vect <- terra::vect(c1[, c(2, 1)], geom = c("lon", "lat"))
    er1 <- terra::rasterize(c1_vect, REF_RAST, fun = sum)
    er1[is.na(er1)] <- 0

    # Calculate ratio
    err <- (er1 / er0)
    # Handle the case of 1/0
    err[is.infinite(err)] <- 1
    # Rest should be NA
    err[err == 0] <- NA
    vr <- terra::values(err)

    # dfr <- terra::as.data.frame(err, xy = TRUE)
    p <- terra::as.points(err)
    data <- cbind(terra::crds(p), terra::values(p))

    return(data)
}

# simulation functions -----

# functions for simulation framework

#' Convert from a named list to a matrix or vector
#' This is just fixing a quirk introduced by a yaml config file
#'
#' @param obj the config object containing specifically phylo config params
#' @param vec logical - if TRUE, converts mean roote vector, if FALSE, converts
#' covar trait matrix
#'
#' @return a matrix of config list elements
#' @export


yamlConvert <- function(obj = config$phylo, vec) {
    n <- obj$ntraits
    if (vec) m <- as.numeric(obj$mean_root_vector)

    if (!vec) {
        nums <- as.numeric(obj$covar_traits)
        m <- matrix(nums, n, n)
    }

    return(m)
}

#' Simulate a tree
#'
#' @param tree_time is amount of time the phylogeny simulation should run
#' around 8 returns 9 tips and 8 interal nodes (for seed 12)
#' birth and death are fixed here
#'
#' @return a tree that was evolved for the specified time
#' @export

fixTree <- function(tree_time) {
    # set.seed(12)
    tr <- ape::rbdtree(0.1, 0, Tmax = tree_time)
    tr$edge.length <- tr$edge.length / max(phytools::nodeHeights(tr)[, 2]) * 1
    return(tr)
}

#' Simulate a trend with phylogeny and traits
#' Simulating an ultrametric tree for the moment
#' Simulating traits by multivariate ornstein-uhlenbeck model
#' Rescaling the tree so that `root depth` is 1, everything else scaled relative
#'
#' @param mean_root_vector mean of traits
#' @param covar_traits covariance among traits
#' @param sigma2 variation in trait
#' @param alpha rubberband parameter
#'
#' @return list with simulated tree and traits
#' @export

traitTrend <- function(tr, ntraits,
                       mean_root_vector,
                       covar_traits,
                       sigma2,
                       alpha) {
    tr_table <- as.data.frame(tr$edge)
    colnames(tr_table) <- c("node_from", "node_to")
    tr_table$length <- tr$edge.length

    trait_matrix <- matrix(0, ntraits, nrow = nrow(tr_table))
    root_num <- tr_table[1, 1]
    root_traits <- MASS::mvrnorm(1, mean_root_vector, covar_traits)

    for (i in seq_len(nrow(tr_table))) {
        if (tr_table[i, 1] == root_num) {
            next_traits <- alpha * root_traits +
                MASS::mvrnorm(1,
                    mu = rep(0, ntraits),
                    Sigma = diag(sigma2 * tr_table[i, "length"], nrow = ntraits)
                )
            trait_matrix[i, ] <- next_traits
        } else {
            old_state_num <- which(tr_table[i, 1] == tr_table[, 2])
            old_state <- trait_matrix[old_state_num, ]
            next_traits <- alpha * old_state +
                MASS::mvrnorm(1,
                    mu = rep(0, ntraits),
                    Sigma = diag(sigma2 * tr_table[i, "length"], nrow = ntraits)
                )
            trait_matrix[i, ] <- next_traits
        }
    }

    trait_matrix <- rbind(root_traits, trait_matrix)
    colnames(trait_matrix) <- paste0("trait-", seq_len(ncol(trait_matrix)))
    node_values <- c(root_num, tr_table$node_to)
    trait_df <- as.data.frame(trait_matrix)
    trait_df$node <- node_values

    return(list(trait_df = trait_df, tr = tr))
}


#' Small function to find distinguish node/tip numbers and split a df
#'
#' @param table values of each internal node and tip. This table MUST have a
#' column named `node` that has node/tip numbers
#' @param tree corresponding phylogenetic tree
#'
#' @return a tip table and node table
#' @export

splitByNodesTips <- function(table, tree) {
    Ntips <- length(tree$tip.label)
    Nnodes <- tree$Nnode

    tip_index <- which(table$node %in% 1:Ntips)
    node_index <- which(!(seq_along(nrow(table)) %in% tip_index))

    tip_table <- table[tip_index, ]
    node_table <- table[node_index, ]

    return(list(tip_table = tip_table, node_table = node_table))
}



# #' Plot ancestral and tip states as branch colours
# #' This is a dumb function because trait has to == be 1/2
# #'
# #' @param tree phylogenetic tree
# #' @param tableOfStuff table with trait values and node numbers
# #' @param trait 1 or 2
# #'
# #' @importFrom ggtree ggtree
# #' @return a plot
# #' @export

# plotTreeBranchColors <- function(tree, tableOfStuff, trait,
# ylow = 0, yhigh = 5) {

#   library(wesanderson)

#   tree2 <- treeio::full_join(tree, tableOfStuff, by = "node")
#   #trait_index <- which(colnames(tree2@data) == traitName)
#   #print(trait_index)

#   if (trait == 1) {
#     pal <- wes_palette("Zissou1", 10, type = "continuous")
#     g <- ggtree::ggtree(tree2, aes(color = tree2@data$`trait-1`),
#     continous = TRUE, size = 2) +
#       #scale_color_distiller(palette = "Spectral", direction = 1) +
#       scale_color_gradientn(colours = pal, limits = c(-8, 9.5)) +
#       geom_tiplab(hjust = -.2)

#   }

#   if(trait == 2) {
#     myPalette <- colorRampPalette(rev(RColorBrewer::brewer.pal(11,
#     "Spectral")))
#     #sc <- scale_colour_gradientn(colours = myPalette(100),
#     #limits=c(ylow, yhigh))
#     pal <- wes_palette("Zissou1", 10, type = "continuous")


#     g <- ggtree::ggtree(tree2, aes(color = tree2@data$`trait-2`),
#      continous = TRUE, size = 2) +
#       scale_color_gradientn(colours = pal, limits = c(-8, 9.5)) +
#       #scale_color_distiller(palette = "Spectral", direction = 1) +
#       geom_tiplab(hjust = -.2) #+ sc
#   }

#   return(g)

# }


#' Simulate data from traits and trees
#'
#' @param n number of samples
#' @param continuous continuous Y data or p/a
#'
#' @return list with simulated tree, traits, XData and Ydata
#' @export

simulateAbundanceData1 <- function(TREE, simulation_params) {
    # set.seed(11)
    # n <- simulation_params$n
    tree_time <- simulation_params$tree_time
    ntraits <- simulation_params$ntraits
    mean_root_vector <- simulation_params$mean_root_vector
    covar_traits <- simulation_params$covar_traits
    sigma2 <- simulation_params$sigma2$value
    alpha <- simulation_params$alpha$value

    trait_tree <- traitTrend(
        tr = TREE, ntraits = ntraits,
        mean_root_vector = mean_root_vector, # c(1,1),
        covar_traits = covar_traits,
        # matrix(c(0.5, 1, 1, 2),2),
        sigma2 = sigma2, alpha = alpha
    )

    tr <- trait_tree$tr
    tr$tip.label <- LETTERS[seq_len(length(tr$tip.label))]
    trait_tbl <- trait_tree$trait_df
    # plotTreeBranchColors(tr, trait_tbl, trait = 2,
    # ylow = min(trait_tbl$`trait-2`),
    # yhigh = max(trait_tbl$`trait-2`))

    ns <- length(tr$tip.label)
    split_tbl <- splitByNodesTips(trait_tbl, tr)
    tip_traits <- split_tbl$tip_table
    rownames(tip_traits) <- tr$tip.label

    return(list(tip_traits = tip_traits, trait_tbl = trait_tbl, tr = tr))
}

#' A second way to simulate abundance
#' Simulating data from trees and traits
#'
#' @param n number of samples
#' @param tip_traits a list of tip trait values
#' @param tr phylogenetic tree
#' @param continuous (logical) whether responses should be continuous or
#' presence/absences
#'
#' @return a list with responses, simulated X data, and true betas
#'
#' @export

simulateAbundanceData2 <- function(n, tip_traits, tr,
                                   continuous = TRUE) {
    # simulating environmental and species data
    n <- n
    habitat <- factor(sample(x = c("forest", "open"), size = n, replace = TRUE))
    climate <- rnorm(n)

    nc <- 4
    ns <- length(tr$tip.label)
    mu <- matrix(0, nrow = nc, ncol = ns)

    # expected niche of each species related to the "covariate" intercept
    mu[1, ] <- -tip_traits$`trait-1`^2 / 4 - tip_traits$`trait-2`
    # expected niche of each species related to the covariate habitat
    mu[2, ] <- 2 * tip_traits$`trait-2`
    # expected niche of each species related to the covariate climate
    mu[3, ] <- tip_traits$`trait-1` / 2
    # expected niche of each species related to the covariate climate*climate
    mu[4, ] <- -1 / 4
    beta <- mu + 0.25 * matrix(rnorm(n = ns * nc), ncol = ns)
    # X = cbind(rep(1, n),habitat_var, climate_var)
    X <- cbind(
        rep(1, n), as.numeric(habitat == "forest"),
        climate, climate * climate
    )
    L <- X %*% beta

    Y <- L + MASS::mvrnorm(n = n, mu = rep(0, ns), Sigma = diag(ns))
    Y2 <- 1 * (L + MASS::mvrnorm(n = n, mu = rep(0, ns), Sigma = diag(ns)) > 0)
    colnames(Y2) <- colnames(Y) <- tr$tip.label

    # putting data together

    XData <- data.frame(climate = climate, habitat = habitat)

    if (continuous == TRUE) Y_return <- Y
    if (!continuous) Y_return <- Y2

    return(list(Y = Y_return, XData = XData, true_betas = beta))
}

#' simulate abundance data somewhat naively - this function just implements
#' a simpler (less noise, more connected) simulation framework for generating
#' abundance data from
#'
#' @param n number of occurrence records
#' @param tip_traits matrix of trait values for tips of the phylogenetic tree
#' @param tr phylogenetic tree
#' @param continuous logical - whether or not to simualte occurrence
#' records (TRUE) or presence-absence (FALSE)
#' @param addNoise2Beta logical - whether or not to jitter betas (TRUE)
#' with white gaussian noise or not
#' @param ind logical - simulate independent betas (TRUE) or not (FALSE)
#' @param data_def logical - whether or not to simulate data-deficiency or not
#' @param def_n number of data deficient records to simulate
#' @param split logical - TRUE if data-deficient species should be split of tree
#' or not
#' @param focalsp character - which species to make data-deficient
#'
#' @return list of split responses, simulated covariate data, split true
#' betas, split tree, full responses, full tree, full true betas, true
#' (complete) Y and new Xdata for prediction
#'
#' @export

simulateAbundanceDataNaive <- function(n,
                                       tip_traits,
                                       tr,
                                       continuous = TRUE,
                                       addNoise2Beta = FALSE,
                                       ind,
                                       data_def,
                                       def_n,
                                       split,
                                       focalsp = "t1") {
    # simulating environmental and species data
    habitat <- rnorm(n)
    climate <- rnorm(n)

    nc <- 1
    ns <- length(tr$tip.label)
    mu <- matrix(0, nrow = nc + 1, ncol = ns)

    # expected niche of each species related to the "covariate" intercept
    mu[1, ] <- 0.25
    # expected niche of each species related to the covariate habitat
    mu[2, ] <- tip_traits$`trait-2` * 5
    # expected niche of each species related to the covariate climate
    # mu[3, ] <- tip_traits$`trait-1`
    # expected niche of each species related to the covariate climate*climate
    # mu[4, ] <- -1/4
    # beta = mu #+ 0.25*matrix(rnorm(n = ns*nc), ncol = ns)
    # X = cbind(rep(1, n),habitat_var, climate_var)
    L <- ape::vcv(tr, corr = TRUE)
    beta <- apply(mu, 1, function(x) {
        MASS::mvrnorm(n = 1, mu = x, Sigma = L)
    })
    if (ind) {
        I <- matrix(0, ns, ns)
        diag(I) <- 1
        beta <- apply(mu, 1, function(x) {
            MASS::mvrnorm(n = 1, mu = x, Sigma = I)
        })
    }

    if (addNoise2Beta) beta + 0.25 * matrix(rnorm(n = ns * nc), ncol = nc)
    X <- cbind(rep(1, n), habitat)
    M <- X %*% t(beta)
    Y <- M + MASS::mvrnorm(n = n, mu = rep(0, ns), Sigma = diag(ns))
    Y2 <- 1 * (M + MASS::mvrnorm(n = n, mu = rep(0, ns), Sigma = diag(ns)) > 0)
    colnames(Y2) <- colnames(Y) <- tr$tip.label

    # putting data together
    XData <- data.frame(habitat = habitat)
    intercept <- rep(1, nrow(XData))
    XData <- cbind(intercept, XData)
    habitat2 <- rnorm(n)
    climate2 <- rnorm(n)
    XData_new <- data.frame(
        intercept = intercept,
        habitat = habitat2
    )

    if (continuous == TRUE) Y_return <- Y
    if (!continuous) Y_return <- Y2

    if (data_def) {
        focal_sp <- focalsp
        index_sp <- which(colnames(Y_return) == focal_sp)
        trueY <- Y_return[, focal_sp]
        Y_return[sample(seq_len(nrow(Y_return)), nrow(Y_return) - def_n,
            replace = FALSE
        ), index_sp] <- 0
    } else {
        trueY <- NULL
    }

    if (split) {
        # drop tip from tr where tip name = t1
        # remove betas and X and Y data for t1
        # save the true tree, and original data deficient X and Y separately

        tr_split <- ape::drop.tip(tr, focal_sp)
        index_sp <- which(colnames(Y_return) == focal_sp)
        Y_return_split <- Y_return[, -c(index_sp)]
        # Y_return_split <- subset(Y_return, select = -c(t1))
        index <- which(rownames(beta) == focal_sp)
        beta_split <- beta[-c(index), ]
    } else {
        tr_split <- tr
        Y_return_split <- Y_return
        beta_split <- beta
    }

    # Everything with "_split" is modified in the split if/else.
    # If split == FALSE,
    # no modification. If TRUE, the focal species has been dropped from
    # the true betas,
    # the Ys and the tree. The full true betas, Ys and the tree are then
    # stored in "_full"
    # trueY is the Y of the data deficient species before it was replaced with 0s

    return(list(
        Y = Y_return_split, XData = XData, true_betas = beta_split,
        tr = tr_split, Y_full = Y_return, true_betas_full = beta,
        tr_full = tr, trueY = trueY, XData_new = XData_new
    ))
}

# TODO: Rewrite this function
#' Check if tree has number of tips in the right range
#' If check fails, function resimulates tree
#'
#' @param a phylo tree
#' @param phylo a list with elements "atleastN" (min number of tips) and
#' "notMoreThan" (max number of tips)
#'
#' @return a tree with tips within the desired range
#' @export

checkTree <- function(TREE1 = tree, PHYLO = phylo) {
    # CASE 1 : AtleastN

    if (!is.null(phylo$atleastN) && is.null(phylo$notMoreThan)) {
        startTime <- Sys.time()
        while (length(tree$tip.label) < phylo$atleastN) {
            tree <- fixTree(tree_time = phylo$tree_time)
            # print("yep")
            endTime <- Sys.time()
            if (as.numeric(endTime - startTime) > 120) {
                print("Warning: you might want to consider changing tree_time instead")
                break
            }
        }
    }

    # CASE 2 : NotMoreThan
    if (is.null(phylo$atleastN) && !is.null(phylo$notMoreThan)) {
        startTime <- Sys.time()
        while (length(tree$tip.label) > phylo$notMoreThan) {
            tree <- fixTree(tree_time = phylo$tree_time)
            # print("yep")
            endTime <- Sys.time()
            if (as.numeric(endTime - startTime) > 120) {
                print("Warning: you might want to consider changing tree_time instead")
                break
            }
        }
    }

    # CASE 3 : BOTH
    if (!is.null(phylo$atleastN) && !is.null(phylo$notMoreThan)) {
        # TODO: Figure this out later
    }

    return(tree)
}

# TODO: Rewrite this function
#' A function to make parameter table for trees and tips
#'
#' @param CONFIG a list with arguments for data_simulation
#'
#' @importFrom emdbook lseq
#' @return a matrix
#' @export

makeParameterTbl <- function(CONFIG = config) {
    PARENT <- CONFIG$data_sim$var$parent
    FT <- CONFIG$data_sim$var$feature
    w <- CONFIG[[PARENT]]
    ind <- which(names(w) %in% FT)

    # Check if we need to interpolate and then interpolate
    keep <- list()
    for (o in seq_len(length(ind))) {
        mm <- w[[ind[o]]]
        if (mm$interpol && length(mm$value) == 2) {
            s <- seq(mm$value[[1]], mm$value[[2]], mm$interpol_step)
            if (mm$log) s <- emdbook::lseq(mm$value[[1]], mm$value[[2]], 10)
        } else {
            s <- as.numeric(c(mm$value))
        }
        keep[[o]] <- s
    }

    d1 <- expand.grid(keep)
    colnames(d1) <- names(w[ind])
    return(d1)
}

#' Clean environment matrix for a given set of variables
#' This function removes any rows that aren't complete
#'
#' @param env matrix
#' @param vars2 selected variables
#'
#' @return a cleaned matrix with no rows with missing data and only selected
#' variables
#'
#' @export

cleanEnv <- function(env1, vars2) {
    w <- which(complete.cases(env1) == FALSE)
    env1 <- env1[-w, ]
    if (any(grepl("\\.tif", vars2))) {
        ind <- grep("\\.tif", vars2)
        for (i in ind) vars2[i] <- gsub("\\.tif", "", vars2[i])
    }
    env1 <- as.matrix(env1[, c(vars2)])
    return(env1)
}

#' Make X matrix
#' This function takes in a matrix and scales it and adds and intercept column
#'
#' @param env1 matrix
#'
#' @return scaled matrix with an intercept column
#'
#' @export

makeXmat <- function(env1) {
    # Make X matrix
    x <- scale(env1)
    # int <- rep(1, nrow(x))
    # e <- rnorm(nrow(x))
    xmat <- matrix(cbind(x), ncol = ncol(x))
    colnames(xmat) <- c(colnames(x))
    return(xmat)
}

#' Give absence - presences indices
#' For a given focal species and a observation data-frame/matrix, it generates
#' two sets of indices - one for the rows where all the focal species are
#' absent and a second for rows where one or more of the focal species are
#' present
#' @param focal_sp a character vector of focal species
#' @param Y large (pool of) observation data frame or matrix
#'
#' @return A list of two elements, one with absence indices and second of
#' presence indicies
#'
#' @export

give01Index <- function(focal_sp, Y) {
    ff <- focal_sp
    ff_len <- length(ff)
    ltab0 <- Y[, ff] == 0
    ltab1 <- Y[, ff] == 1

    if (length(ff) > 1) {
        fsp_abs <- which(rowSums(ltab0) == ff_len)
        fsp_pres <- which(rowSums(ltab1) >= 1)
    } else {
        fsp_abs <- which(ltab0 == ff_len)
        fsp_pres <- which(ltab1 == 1)
    }

    return(list("index0" = fsp_abs, "index1" = fsp_pres))
}

#' create a dataset with specific absence and presence specs
#'
#' For a given observation pool matrix, list of indices, number of samples and
#' number of data-deficient species samples, create a data with dd specs
#'
#' @param Y large (pool of) observation matrix
#' @param ind list of presence/absence indices named "index0" and "index1"
#' @param nsamples number of total samples desired
#' @param def_nsamples number of total data-deficienct samples desired
#'
#' @return final observation matrix with nrow = nsamples and focal species
#' will have presences = def_nsamples
#'
#' @export

create01Dataset <- function(Y, ind, nsamples, def_nsamples) {
    fsp_abs <- ind$index0
    fsp_pres <- ind$index1

    # Absence data
    sub0 <- Y[fsp_abs, ]
    sub0 <- sub0[order(sub0$tot, decreasing = TRUE), ]
    grab0 <- nsamples - def_nsamples
    if (grab0 > nrow(sub0)) {
        warning("Not enough absence samples")
    }
    sub0 <- sub0[seq_len(grab0), ]

    # Presence data
    sub1 <- Y[fsp_pres, ]
    sub1 <- sub1[order(sub1$tot, decreasing = TRUE), ]
    grab1 <- def_nsamples
    if (grab1 > nrow(sub1)) {
        warning("Not enough presence samples")
    }
    sub1 <- sub1[seq_len(grab1), ]
    Y_fin <- rbind(sub0, sub1)
    return(Y_fin)
}

#' Print details
#' Generic function but specialized here to print the dimensions and column
#' sums of observed test and train dataframes
#'
#' @param digits digits
#' @param width width
#' @param ... optional arguments
#'
#' @return printed summary
#'
#' @export

printDetails <- function(..., digits = getOption("digits"),
                         width = getOption("width")) {
    op <- options(digits = digits, width = width)
    on.exit(options(op))
    call <- match.call()
    call[c("digits", "width")] <- NULL

    f1 <- function(x, name) {
        paste0(
            "For ", name, " specs are as follows: Dimensions are ",
            paste0(format(dim(x)), collapse = ", "),
            " and number of presence points per species are ",
            paste0(format(colSums(x)), collapse = ", "), "."
        )
    }
    l <- Map(f1,
        x = list(...), name = lapply(call[-1], deparse),
        USE.NAMES = FALSE
    )
    s <- do.call(paste, c(l, list(sep = "\n\n")))
    writeLines(strwrap(s))
}

#' Create either a systematic sample/biased grid
#'
#' Based on user input create either a grid of systematically sampled
#' coordinates or biased (user-input) coordinates
#'
#' @param user_grid Should a user (biased) grid be created (TRUE) or a
#' systematically sampled grid (FALSE)
#' @param nsamples Number of samples desired in resulting dataset. This function
#' takes it in and multiples by 20 to create a 'sufficiently' large pool to
#' create data-deficient matricies easily downstream
#' @param path relative path to the user grid
#'
#' @return a large grid of coordinates
#'
#' @export

createGrid <- function(user_grid,
                       nsamples = 100,
                       path,
                       xmin = -10,
                       xmax = 10,
                       ymin = -10,
                       ymax = 10) {
    #### Sample a systematic grid for covariates ####
    pool <- nsamples * 50
    if (!user_grid) {
        pool2 <- ceiling(sqrt(pool))
        x_seq <- seq(xmin, xmax, length.out = pool2)
        y_seq <- seq(ymin, ymax, length.out = pool2)
        sys_cood <- expand.grid(x = x_seq, y = y_seq)
        grid <- sys_cood
    } else {
        # Use user coordinates for sampling grid
        # Load corresponding raw data
        load(file = path)
        # loads object cood
        user_cood <- cood[, c("x", "y")]
        # this is so dumb buy have to reorder x and y because apparently
        # terra isn't smart enough to get it
        rm(cood)
        if (pool > nrow(user_cood)) error("Not enough user input coordinates")
        grid <- user_cood[sample(nrow(user_cood), size = pool), ]
    }

    return(grid)
}

#' wrangle environmental data
#' Read in environmental data based on specified variables and crop and clean
#'
#' @param vars1 environmental variables
#' @param extent vector of xmin, xmax, ymin, ymax - not used
#' @return cleaned and cropped and wrangled environmental data matrix
#'
#' @importFrom terra rast
#' @importFrom terra crop
#' @importFrom terra extract
#'
#' @export

wrangleEnv <- function(vars1 = c("CHELSA_bio_1.tif", "CHELSA_bio_13.tif"),
                       extent = c(-81, -60, -40, 0)) {
    ##### Read in env covariates ####
    # just temp and mean annual precip and elev - Bio1 and Bio13
    epath <- "~/env"

    env_files <- file.path(epath, vars1)
    ras <- terra::rast(env_files)

    # Grab data points
    env <- terra::extract(ras, grid)
    env <- cleanEnv(env, vars2 = vars1)

    return(env)
}

#' Simulate Abundance Data in a relatistic way
#'
#' @param grid grid of coordinates
#' @param data_deficiency make something data-deficient (yes/TRUE)
#' @param focal_sp character vector of focal species
#' @param nsamples number of total desired samples (note these will be split
#' 70/30 as training/test datasets)
#' @param binary continuous data (FALSE) or presence-absence (TRUE)
#' @param print print out details (yes/TRUE)
#' @param tip_traits result from simulate abundance 1 function
#' @param tree corresponding phylogenetic tree
#' @param OUbetas (logical) TRUE if OU function should be used to simulate
#' betas, FALSE if betas should be mvrnorm using traitTrend function
#' @param params phylogenetic parameters
#' @param PRACTICE_WHEELS true or false - true for GP covar simulation
#' @param vars variables to pull environmental data from for the simulation
#' framework
#'
#' @return a list of data
#'
#' @export

simulateAbundanceData_irl <- function(grid,
                                      data_deficiency = TRUE,
                                      focal_sp = "A",
                                      nsamples = 100,
                                      def_nsamples = 30,
                                      binary = TRUE,
                                      print = TRUE,
                                      tip_traits,
                                      tree,
                                      OUbetas = TRUE,
                                      params = NULL,
                                      PRACTICE_WHEELS,
                                      vars = c("CHELSA_bio_1.tif", "CHELSA_bio_13.tif")) { # nolint
    # TODO: i) Make systematic grid versus user grid an option
    # ii) data-deficiency options
    # iii) class imbalance options

    env <- wrangleEnv(vars1 = vars)

    # Blank names
    sp <- tree$tip.label

    #### Simulate responses ####
    # Simulate with real env data
    nc <- ncol(tip_traits) - 1
    ns <- length(tree$tip.label)

    if (OUbetas) {
        if (is.null(params)) error("Please specify simulation parameters")
        beta <- OUCovar(TREE = tree, PARAMS = params)
    } else {
        if (PRACTICE_WHEELS) {
            if (is.null(params)) error("Please specify simulation parameters")
            beta <- practice_wheels_covar(TREE = tree, PARAMS = params)
        } else {
            # +1 for intercept - intercept removed
            mu <- matrix(0, nrow = (nc), ncol = ns)

            # Expected value
            # mu[1, ] <- 0.2 # Intercept
            mu[1, ] <- tip_traits$`trait-1`
            mu[2, ] <- tip_traits$`trait-2`

            # Note vcv assumes BM covariance structure
            # For OU a correlation structure needs to be specified
            # L <- ape::vcv(tree)
            beta <- apply(mu, 1, function(x) {
                MASS::mvrnorm(n = 1, mu = x, Sigma = diag(1, ns, ns))
            })
        }
    }

    # beta <- beta + 3

    # Systematic sampling grid first
    Xmat <- makeXmat(env)
    # Xmat <- Xmat[, -c(which(colnames(Xmat) == "Effort"))]
    M <- Xmat %*% t(beta)
    pool2 <- nrow(Xmat)

    # Independently simulating species observations
    covar <- diag(ns)
    diag(covar) <- 2
    Y <- M + MASS::mvrnorm(n = nrow(M), mu = rep(0, ns), Sigma = covar)

    # Fineaggling
    Y2 <- matrix(as.numeric(Y > 0), ncol = ncol(Y))
    Y2 <- data.frame(Y2)
    Y <- data.frame(Y)
    colnames(Y2) <- colnames(Y) <- sp
    Xmat <- data.frame(Xmat)
    Y2$rwid <- Xmat$rwid <- Y$rwid <- paste0("rw_", seq_len(pool2))

    # Create temporary frame so Y can remain the "groud truth"
    Y_temp <- Y2
    Y_temp$tot <- rowSums(Y_temp[, sp])

    if (data_deficiency) {
        # Data deficiency
        # First let's create various indexes based on categories of Y
        ind <- give01Index(focal_sp, Y_temp)

        # Create a data_deficiency
        Y_fin <- create01Dataset(Y = Y_temp, ind, nsamples, def_nsamples)
        test1 <- which(!complete.cases(Y_fin))

        # Check
        if (identical(length(test1), 0L) == FALSE) {
            warning("There are NAs in the simulated Y data")
        }

        test2 <- colSums(Y_fin[, sp])
        if (any(test2[focal_sp] > def_nsamples)) {
            warning("More than required data-def samples")
        }

        if (any(test2[which(names(test2) %in% focal_sp == FALSE)] < def_nsamples)) {
            warning("Some abundant species are accidentally deficient")
        }
    } else {
        k <- sample(seq_len(nrow(Y_temp)), size = nsamples)
        Y_fin <- Y_temp[k, ]
    }

    # More checks
    X_fin <- Xmat[match(Y_fin$rwid, Xmat$rwid), ]
    Yc_fin <- Y[match(Y_fin$rwid, Y$rwid), ]

    # Another check
    if (any(X_fin$rwid != Y_fin$rwid)) {
        warning("Misalignment in X and Y dataframes")
    }

    if (any(nrow(X_fin) != nrow(Y_fin))) {
        warning("Number of rows in X and Y differ")
    }

    if (nrow(Y_fin) != nsamples) {
        warning("Number of observations != number of samples")
    }

    if (nrow(Y2) != nrow(Xmat)) {
        warning("The full datasets are mismatched")
    }

    # Create train and test datasets
    ntrain <- round(0.7 * nsamples, 0)
    ntest <- nsamples - ntrain
    ind_train <- sample(seq_len(nsamples), ntrain)
    Y_train <- Y_fin[ind_train, ]
    Yc_train <- Yc_fin[ind_train, ]
    Y_test <- Y_fin[-ind_train, ]
    Yc_test <- Yc_fin[-ind_train, ]
    X_train <- X_fin[ind_train, ]
    X_test <- X_fin[-ind_train, ]

    if (binary) {
        Yr <- Y2
        Yr_train <- Y_train
        Yr_test <- Y_test
    } else {
        Yr <- Y
        Yr_train <- Yc_train
        Yr_test <- Yc_test
    }

    if (print) {
        printDetails(Y_train[, sp], Y_test[, sp], digits = 4L, width = 72L)
    }

    e <- list(
        "Ytrue" = Yr, "X" = Xmat, "Ytrain" = Yr_train, "Xtrain" = X_train,
        "true_betas" = beta, "train_index" = ind_train,
        "Ytest" = Yr_test, "Xtest" = X_test
    )
    return(e)
}

#' A function to simulate beta values from an OU process
#' Using the covariance function between two particles under an
#' OU process, the covariance between tips is calcualted and values
#' simulated based on this covariance
#'
#' @param TREE phylogenetic tree
#' @param PARARMS list of phylogenetic parameters that contain elements
#' sigma2$value and alpha$value and ntraits
#'
#' @return a vector of betas the same length as number of tips in the
#' phylogenetic tree
#'
#' @importFrom MASS mvrnorm
#' @importFrom ape cophenetic.phylo
#'
#' @export

OUCovar <- function(TREE, PARAMS) {
    # Collect data and set up empty covariance matrix
    sigma2 <- PARAMS$sigma2$value
    K <- PARAMS$alpha$value
    dist <- ape::cophenetic.phylo(TREE)
    ntips <- length(TREE$tip.label)
    covar <- matrix(0, ncol = ncol(dist), nrow = nrow(dist))
    ntraits <- PARAMS$ntraits

    # OU covariance
    for (i in 1:ntips) {
        for (j in 1:ntips) {
            covar[i, j] <- (sigma2 / (2 * K)) * exp(-K * dist[i, j])
        }
    }

    means <- c(sample(c(3, 4), ntraits, replace = TRUE))

    betas <- matrix(0, ncol = length(means), nrow = ntips)
    for (i in seq_len(length(means))) {
        est <- MASS::mvrnorm(n = 1000, mu = c(rep(means[i], ntips)), covar)
        betas[, i] <- colMeans(est)
    }

    return(betas)
}

#' A function to simulate covariance using a gaussian process
#' Practice wheels because it's easier for the GP model to estimate (in theory)
#' @param TREE phylogenetic tree
#' @param PARARMS list of phylogenetic parameters that contain elements
#' rho$value and alphaw$value and ntraits
#'
#' @return a vector of betas the same length as number of tips in the
#' phylogenetic tree
#'
#' @importFrom MASS mvrnorm
#' @importFrom ape cophenetic.phylo
#'
#' @export

practice_wheels_covar <- function(TREE, PARAMS) {
    dist <- ape::cophenetic.phylo(TREE)
    ntips <- length(TREE$tip.label)
    covar <- matrix(0, ncol = ncol(dist), nrow = nrow(dist))
    ntraits <- PARAMS$ntraits
    rho <- PARAMS$rho$value
    alpha2 <- PARAMS$alpha2$value

    # GP covariance
    for (i in 1:ntips) {
        for (j in 1:ntips) {
            covar[i, j] <- alpha2^2 * exp(-0.5 * (dist[i, j]^2) / rho^2)
        }

        means <- c(1, sample(c(2, 3), ntraits, replace = TRUE)) # 1 for intercept

        betas <- matrix(0, ncol = length(means), nrow = ntips)
        for (i in seq_len(length(means))) {
            est <- MASS::mvrnorm(n = 1000, mu = c(rep(means[i], ntips)), covar)
            betas[, i] <- colMeans(est)
        }
    }
    return(betas)
}


# 010-raw_data_functions -----

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

# organizeCovar <- function(epath, env.files, env.crs, cood) {
#     env.paths <- file.path(epath, env.files)
#     colnames(cood) <- c("y", "x")
#     cood <- cood[, c("x", "y")]
#     extent <- raster::extent(cood)

#     env <- raster::stack()
#     for (env.path in env.paths) {
#         print(paste0("Loading ", env.path, " ..."))
#         env.layer <- raster::raster(env.path)

#         # If needed, reproject the env layer
#         if (!raster::compareCRS(raster::crs(env.layer), env.crs)) {
#             print(paste0(" Reprojecting..."))
#             env.layer <- raster::projectRaster(env.layer, env.crs)
#         }

#         if (!is.null(extent)) {
#             env.crop <- raster::crop(env.layer, extent)
#         }
#         env <- raster::addLayer(env, env.crop)
#     }

#     cood.mat <- do.call(cbind, cood)

#     annotations <- raster::extract(env, cood.mat,
#         method = "bilinear"
#     )
#     return(annotations)
# }

# Function for creating a reference coord set

makeRefCoods <- function(REF_RAST, AREA = area) {
    REF_RAST[] <- 1
    area_vect <- terra::vect(AREA)
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


createEffortRast <- function(SAMPLE_RAST, EC, VAR, START_TIME = 2003,
                             END_TIME = NULL, PATH = NULL) {
    if (!is.null(END_TIME)) {
        EC <- EC[which(EC$year >= START_TIME & EC$year <= END_TIME), ]
    } else {
        EC <- EC[which(EC$year >= START_TIME), ]
    }
    e1_vect <- terra::vect(EC[, c("lat", "lon", VAR)],
        geom = c("lon", "lat")
    )
    e1 <- terra::rasterize(e1_vect, SAMPLE_RAST, field = VAR, fun = sum)

    if (!is.null(PATH)) {
        if (!is.null(END_TIME)) {
            name <- paste0("effort_rast_", VAR, "_e", paste0(START_TIME, "-e", END_TIME), ".tiff")
        } else {
            name <- paste0("effort_rast_", VAR, "_FULL.tiff")
        }
        terra::writeRaster(e1, file = file.path(PATH, name), overwrite = TRUE)
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
        # print(paste0("Loading ", env_path, " ..."))
        env_layer <- terra::rast(env_path)

        # Reproject
        # print(paste0(" Reprojecting..."))
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
        # print(i)
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

#' A map of occurrence records for a selected species
#'
#' Creates a map of Presence/Absence points for a selected species
#'
#' @param YC occurrence dataset (y) with coordinates (c). The coordinates must
#' have columns "x" and "y"
#' @param FSP (f)ocal )sp)ecies i.e. the name of the species of interest
#'
#' @importFrom rnaturalearth ne_countries
#' @importFrom raster extent
#' @import ggplot2
#'
#' @return A ggplot object (map in this case)
#'
#' @export

makeOccurrenceMap <- function(YC, FSP, WORLD = world) {
    if (all(c("lat", "lon") %in% colnames(YC))) {
        ind_lat <- which(colnames(YC) == "lat")
        ind_lon <- which(colnames(YC) == "lon")
        colnames(YC)[c(ind_lon, ind_lat)] <- c("x", "y")
    }
    COOD <- YC[, c("x", "y")]
    ex <- getExtentDf(COOD, NAMES = c("x", "y"))
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
            data = WORLD, fill = NA,
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
    er0 <- terra::rasterize(c0_vect, REF_RAST, fun = sum)
    er0[is.na(er0)] <- 0
    c1 <- YCOOD[ind1, ]
    c1_vect <- terra::vect(c1[, NAMES], geom = NAMES)
    er1 <- terra::rasterize(c1_vect, REF_RAST, fun = sum)
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
        if (length(ind1) == 0) {
            print("No data")
            data[[j]] <- NA
        } else {
            # Returns the ratio of presences to absences at every pixel
            err_df <- gimmeRatios(yfsp, YC_CROP[, c("x", "y")], ER)
            err_df$sp <- SPS[j]
            data[[j]] <- err_df
        }
    }
    return(data)
}

createBgPoints <- function(
    COOD, AREA_VECT, EX, NROW_Y, NCOL_Y, MULTIPLE = 1, BUFFER = 10000,
    EXPAND_FRAC = 0.1) {
    # Convert the input coordinates to a SpatVector with the specified CRS
    cood_vect <- terra::vect(COOD, geom = c("lon", "lat"), crs = env_crs)

    # Calculate expansion in x and y directions based on the extent and expansion fraction
    dx <- (EX["xmax"] - EX["xmin"]) * EXPAND_FRAC
    dy <- (EX["ymax"] - EX["ymin"]) * EXPAND_FRAC

    # Create a widened extent around the original extent
    ex_wide <- terra::ext(
        EX["xmin"] - dx,
        EX["xmax"] + dx,
        EX["ymin"] - dy,
        EX["ymax"] + dy
    )

    # Crop the area vector to the widened extent
    area_vect_sp <- terra::crop(AREA_VECT, ex_wide)

    # Buffer the coordinates by the specified width (in meters)
    cood_vect_buffer <- terra::buffer(cood_vect, width = BUFFER)

    # If the buffered area is larger than the cropped area, expand the area further
    if (terra::ext(area_vect_sp) < terra::ext(cood_vect_buffer)) {
        new_ext <- terra::ext(cood_vect_buffer)
        dx <- (terra::xmax(new_ext) - terra::xmin(new_ext)) * EXPAND_FRAC
        dy <- (terra::ymax(new_ext) - terra::ymin(new_ext)) * EXPAND_FRAC
        new_ext_wide <- terra::ext(
            terra::xmin(new_ext) - dx,
            terra::xmax(new_ext) + dx,
            terra::ymin(new_ext) - dy,
            terra::ymax(new_ext) + dy
        )
        area_vect_sp <- terra::crop(AREA_VECT, new_ext_wide)
    }

    # Remove the buffered area (around the coordinates) from the area vector to avoid overlap
    area_vect_sp <- terra::erase(area_vect_sp, cood_vect_buffer)
    repeat {
        bg_cood_vect <- terra::spatSample(area_vect_sp, (NROW_Y * MULTIPLE), method = "random")

        # bg_cood matrix
        bg_cood <- terra::geom(bg_cood_vect)[, c("x", "y")]
        bg_cood <- data.frame(bg_cood)
        colnames(bg_cood) <- c("lon", "lat")

        x1 <- annotateCoods(epath, env_files, env_crs, COOD = bg_cood)
        xx1 <- annotateCoods(
            effpath,
            eff_files, env_crs, bg_cood
        )
        colnames(xx1) <- c("distance", "duration", "num_observers")
        x1 <- cbind(x1, xx1)
        x1[, "distance"] <- rep(1, nrow(x1)) # will be logged, log(1) = 0
        x1[, "duration"] <- rep(1, nrow(x1))
        x1[, "num_observers"] <- rep(1, nrow(x1))

        y_bg <- matrix(0, nrow = nrow(bg_cood), ncol = NCOL_Y)

        # Only complete cases
        rmv1 <- which(complete.cases(x1) == FALSE)
        if (length(rmv1) > 1) {
            x1 <- x1[-rmv1, ]
            y_bg <- y_bg[-rmv1, ]
            bg_cood <- bg_cood[-rmv1, ]
        }

        # Check if there are at least 2 non-NA rows
        if (nrow(x1) >= 2) {
            break
        }
    }
    return(list("y_bg" = y_bg, "x1" = x1, "bg_cood" = bg_cood))
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

# 999-function_bounce.R -----

# Bounce file for 032-spatial_prediction.R debugging and function writing

# TODO: After testing, all these need to be sent over to phyloGenie
# kmeans for binning continuous probabilities for
kmeansMaxK <- function(DF, init, VAR = "pp") {
    obj <- "empty"
    while (obj == "empty") {
        j <- which(colnames(DF) == VAR)
        ii <- seq(3, init, by = 1)
        out <- list()
        count <- 1
        if (length(unique(DF[, j])) == 1) {
            # print("y")
            ii <- 1
        }
        for (i in ii) {
            out[[count]] <- tryCatch(
                {
                    cluster <- kmeans(DF[, j], i)
                },
                error = function(cond) {
                    return("empty")
                }
            )
            count <- count + 1
            # print(count)
        }
        ind <- which(sapply(out, class) == "kmeans")
        rr <- max(ind)
        if (!is.finite(rr)) {
            warning("Something went wrong in getting index")
            rr <- 1
        }
        obj2 <- out[[rr]]

        if (length(obj2) == 1) {
            obj <- "empty"
        } else {
            obj <- "not_empty"
        }
    }
    return(obj2)
}

# Function for processing all the config parameters we need

grabConfig <- function(ROOT, EXP_NAME) {
    files <- list.files(file.path(ROOT, "res", EXP_NAME))
    obj <- load(file.path(ROOT, "res", EXP_NAME, files[1]))

    config <- all$args
    def_lev <- config$data$def_level
    fsp <- config$data$focal_sp
    rep <- config$data$nrep
    genus <- config$data$cluster_name
    repno <- config$data$repno

    ll <- list(
        "config" = config, "def_lev" = def_lev,
        "fsp" = fsp, "rep" = rep, "genus" = genus, "repno" = repno
    )

    return(ll)
}

# Function to harmonize lat/lon

latLonNameConvertXY <- function(DF) {
    if (any(c("lat", "lon") %in% colnames(DF))) {
        latind <- which(colnames(DF) == "lat")
        if (length(latind) > 0) colnames(DF)[latind] <- "y"
        lonind <- which(colnames(DF) == "lon")
        if (length(lonind) > 0) colnames(DF)[lonind] <- "x"
    }
    return(DF)
}


# Helper functions for the widening extent stuff


chooseMin <- function(NUM, NUMNP) {
    if (any(NUMNP == "neg")) {
        ind <- which(NUMNP == "neg")
        ind2 <- which(abs(NUM[ind]) == max(abs(NUM[ind])))
        cmin <- NUM[ind2]
    } else {
        # Else is everything is positive in which case we just want smallest no
        ind <- which(NUM == min(NUM))
        cmin <- NUM[ind]
    }
    return(cmin)
}

chooseMax <- function(NUM, NUMNP) {
    if (any(NUMNP == "neg")) {
        ind <- which(NUMNP == "neg")
        ind2 <- which(abs(NUM[ind]) == min(abs(NUM[ind])))
        cmax <- NUM[ind2]
    } else {
        # Else is everything is positive in which case we just want biggest no
        ind <- which(NUM == max(NUM))
        cmax <- NUM[ind]
    }
    return(cmax)
}

# Function for widening extent

chooseBigEx <- function(DOMAIN) {
    # The number we want depends on where it falls on the grid and
    # whether it's for the min or max of the extent
    # Essentially it needs to be a more involved script than this

    # Let's figure out what's negative and what's positive here
    tmp <- (DOMAIN == abs(DOMAIN))
    tmp <- ifelse(tmp, "pos", "neg") # for readability

    # Let's tackle xmin first
    # For -ve, we want the bigger negative number
    # For +ve, we want the smaller positive number
    # For mixed, we want the bigger negative number

    cxmin <- chooseMin(DOMAIN[, "xmin"], tmp[, "xmin"])
    # Now we do xmax
    # For -ve, we want the smaller negative number
    # For +ve, we want the bigger positive number
    # For mixed, we want the bigger positive number
    cxmax <- chooseMax(DOMAIN[, "xmax"], tmp[, "xmax"])

    # Let's tackle ymin first
    # For -ve, we want the bigger negative number
    # For +ve, we want the smaller positive number
    # For mixed, we want the bigger negative number

    cymin <- chooseMin(DOMAIN[, "ymin"], tmp[, "ymin"])

    # Now we do ymax
    # For -ve, we want the smaller negative number
    # For +ve, we want the bigger positive number
    # For mixed, we want the bigger positive number

    cymax <- chooseMax(DOMAIN[, "ymax"], tmp[, "ymax"])

    domain_wide <- terra::ext(c(cxmin, cxmax, cymin, cymax))
    return(domain_wide)
}


createDirs <- function(ROOT, EXP_ROOT, GENUS, FSPALL) {
    # Create the file structure here
    # Main cluster files
    main_file <- file.path(ROOT, "analysis", EXP_ROOT, GENUS)
    if (!file.exists(main_file)) dir.create(main_file)

    sub_files <- paste0(main_file, c(
        "/spatial_pred/",
        "/auc/", "/sens/", "/tss/", "/bias/", "/beta/",
        "/cond_pred/", "/figures/"
    ))

    # Cluster file sub dirs
    if (any(dir.exists(sub_files) == FALSE)) {
        ind <- which(!dir.exists(sub_files))
        for (p in ind) dir.create(sub_files[p])
    }

    # Per focal species
    sp_files <- file.path(main_file, FSPALL)

    if (any(dir.exists(sp_files) == FALSE)) {
        ind <- which(!dir.exists(sp_files))
        for (p in ind) dir.create(sp_files[p])
    }

    # Sub dirs per focal species

    sub_sp_files <- paste0(sp_files, rep(
        c(
            "/spatial_pred/",
            "/auc/", "/sens/", "/tss/", "/bias/", "/beta/",
            "/cond_pred/", "/figures/"
        ),
        each = length(sp_files)
    ))

    # Cluster file sub dirs
    if (any(dir.exists(sub_sp_files) == FALSE)) {
        ind <- which(!dir.exists(sub_sp_files))
        for (p in ind) dir.create(sub_sp_files[p])
    }
}

multiEnvIntersect <- function(L) {
    ids <- sapply(L, function(x) {
        return(x[, "ID"])
    })

    in0 <- intersect(ids[[1]], ids[[2]])
    for (i in 3:length(ids)) in0 <- intersect(in0, ids[[i]])

    l2 <- lapply(L, function(x) {
        ids <- which(x[, "ID"] %in% in0)
        x <- x[ids, ]
        return(x)
    })

    return(l2)
}
# Using this space right here to figure out the terra::rast issue

# m <- matrix(1:25, nrow = 5, ncol = 5)
# rastm <- terra::rast(m)

# d <- as.data.frame(m)
# rastd <- terra::rast(d, type = "xyz")

# tmp <- pts_error[1:25, ]
# tmp[1, "pp"] <- 0.5
# tmpr <- terra::rast(tmp, type = "xyz")
# tmpv <- terra::vect(tmp, geom = c("x", "y"))
# tmpr <- terra::rast(tmpv,
#     crs = "+proj=longlat +datum=WGS84", resolution = 0.04166667
# )
# # Resolution specs are from env files processed in prediction_matricies
# z <- terra::rasterize(tmpv, tmpr, "pp")
# df <- terra::as.data.frame(z, xy = TRUE)

# Define the function to create and attach matrices
create_and_attach_matrices <- function(
    tags, corrected1,
    matrix_names, env = .GlobalEnv) {
    for (name in matrix_names) {
        mat <- matrix(0,
            ncol = length(tags), nrow = length(corrected1),
            dimnames = list(corrected1, tags)
        )
        assign(name, mat, envir = env)
    }
}

# Boyce index
calculate_boyce <- function(PRED_TEST) {
    pred_test <- PRED_TEST
    pred_boyce <- ROCR::prediction(pred_test[, "test"], pred_test[, "true"])
    pred_boyce <- ROCR::performance(pred_boyce,
        measure = "tpr", x.measure = "fpr"
    )
    boyce_index <- unlist(lapply(seq_along(pred_boyce@y.values), function(i) {
        mean(pred_boyce@y.values[[i]])
    }))

    return(boyce_index)
}


# Function to calculate True Skill Statistic (TSS)
calculate_tss <- function(predicted, observed) {
    # Ensure binary classification (0 and 1)
    if (!all(predicted %in% c(0, 1)) | !all(observed %in% c(0, 1))) {
        stop("Predicted and observed values must be binary (0 or 1).")
    }

    # Check if the input vectors have more than one value
    if (length(predicted) <= 2 || length(observed) <= 2) {
        return(NA)
    } else {
        # Create a confusion matrix
        confusion_matrix <- table(predicted, observed)

        # Extract values from confusion matrix
        tp <- confusion_matrix[2, 2] # True Positives
        fn <- confusion_matrix[2, 1] # False Negatives
        fp <- confusion_matrix[1, 2] # False Positives
        tn <- confusion_matrix[1, 1] # True Negatives

        # Calculate Sensitivity (True Positive Rate)
        sensitivity <- tp / (tp + fn)

        # Calculate Specificity (True Negative Rate)
        specificity <- tn / (tn + fp)

        # Calculate TSS
        tss <- sensitivity + specificity - 1

        return(tss)
    }
}

process_x_matrix <- function(x_matrix) {
    # Step 1: Add column names
    colnames(x_matrix) <- c(
        "meanTemp", "tempSeasonality", "precipQuart",
        "precipSeasonality", "cloudCover", "Annual_EVI", "TRI", "elevation",
        "distance", "duration", "num_observers"
    )

    # Step 2: Add new predictors
    x_matrix$meanTemp2 <- (x_matrix$meanTemp)^2
    x_matrix$precipQuart2 <- (x_matrix$precipQuart)^2

    # Step 3: Scale the data
    # x_matrix <- scale(x_matrix)

    # Step 4: Convert to data frame and add Intercept
    x_matrix <- as.data.frame(x_matrix)
    x_matrix$Intercept <- rep(1, nrow(x_matrix))

    # Step 5: Convert back to matrix
    x_matrix <- as.matrix(x_matrix)

    # Step 6: Reorder columns
    order <- c(
        "Intercept", "meanTemp", "meanTemp2", "tempSeasonality",
        "precipQuart", "precipQuart2", "precipSeasonality", "elevation",
        "Annual_EVI", "TRI", "cloudCover", "duration", "distance",
        "num_observers"
    )
    x_matrix <- x_matrix[, order]

    # Return the processed data

    return(x_matrix)
}


# A couple of functions
simOU_new <- function(
    TREE,
    NTRAITS, SIGMA2, ALPHA) {
    if (ALPHA == 0) {
        bb <- phytools::fastBM(TREE, sig2 = SIGMA2, a0 = 1, nsim = NTRAITS)
    } else {
        bb <- phytools::fastBM(TREE,
            alpha = ALPHA, theta = 1,
            sigma2 = SIGMA2, nsim = NTRAITS
        )
    }
    colnames(bb) <- paste0("trait-", seq_len(ncol(bb)))
    return(bb)
}

# simulateAbundanceData_new <- function(grid,
#                                       nsamples = 100,
#                                       binary = TRUE,
#                                       print = TRUE,
#                                       tree,
#                                       ntraits,
#                                       sigma2,
#                                       alpha,
#                                       epath,
#                                       vars = c("CHELSA_bio_1.tif", "CHELSA_bio_13.tif")) { # nolint

#     # Get environmental data
#     # env <- wrangleEnv(vars1 = vars, epath1 = epath, grid1 = grid)
#     env <- wrangleEnv(vars1 = vars)

#     # Blank names
#     sp <- tree$tip.label

#     #### Simulate responses ####

#     ns <- length(sp)

#     # +1 for intercept - intercept removed
#     # mu <- matrix(0, nrow = (nc), ncol = ns)

#     # Expected value
#     # mu[1, ] <- 0.2 # Intercept
#     # mu[1, ] <- tip_traits$`trait-1`
#     # mu[2, ] <- tip_traits$`trait-2`

#     # Note vcv assumes BM covariance structure
#     # For OU a correlation structure needs to be specified
#     # L <- ape::vcv(tree)
#     # beta <- apply(mu, 1, function(x) {
#     #   MASS::mvrnorm(n = 1, mu = x, Sigma = diag(1, ns, ns))
#     # # })

#     beta <- simOU_new(
#         tree,
#         ntraits, sigma2, alpha
#     )
#     # beta <- cbind(rep(1, ns), rep(1, ns))

#     # beta <- beta + 3

#     # Systematic sampling grid first
#     Xmat <- makeXmat(env)

#     # Assuming Xmat is of dimensions [N, K] and t(beta) is [K, J]
#     # The resulting M will have dimensions [N, J]
#     M <- Xmat %*% t(beta)
#     pool2 <- nrow(Xmat)

#     # Set up covariance matrix
#     ns <- ncol(M)
#     covar <- diag(ns)
#     diag(covar) <- 2
#     Y <- M + MASS::mvrnorm(
#         n = nrow(M),
#         mu = rep(0, ns), Sigma = covar
#     )

#     # Define bounds for truncation
#     # upper_bound <- rep(Inf, ns) # No upper truncation
#     # lower_bound <- rep(0, ns) # Lower truncation at 0

#     # Sample truncated multivariate normal for each row in M
#     # Y <- matrix(NA, nrow = nrow(M), ncol = ns)
#     # for (i in seq_len(nrow(M))) {
#     #     Y[i, ] <- tmvtnorm::rtmvnorm(
#     #         n = 1,
#     #         mean = M[i, ],
#     #         sigma = covar,
#     #         lower = lower_bound,
#     #         upper = upper_bound,
#     #         algorithm = "gibbs"
#     #     )
#     #     # Y[i, ] <- MASS::mvrnorm(
#     #     #     n = 1,
#     #     #     mu = M[i, ],
#     #     #     Sigma = covar
#     #     # )
#     # }

#     # Fineaggling
#     Y2 <- matrix(as.numeric(Y > 0), ncol = ncol(Y))
#     Y2 <- data.frame(Y2)
#     Y <- data.frame(Y)
#     colnames(Y2) <- colnames(Y) <- sp
#     Xmat <- data.frame(Xmat)
#     Y2$rwid <- Xmat$rwid <- Y$rwid <- paste0("rw_", seq_len(pool2))

#     # Create temporary frame so Y can remain the "groud truth"
#     Y_temp <- Y2
#     Y_temp$tot <- rowSums(Y_temp[, sp])

#     k <- sample(seq_len(nrow(Y_temp)), size = nsamples)
#     Y_fin <- Y_temp[k, ]

#     # More checks
#     X_fin <- Xmat[match(Y_fin$rwid, Xmat$rwid), ]
#     Yc_fin <- Y[match(Y_fin$rwid, Y$rwid), ]

#     # Another check
#     if (any(X_fin$rwid != Y_fin$rwid)) {
#         warning("Misalignment in X and Y dataframes")
#     }

#     if (any(nrow(X_fin) != nrow(Y_fin))) {
#         warning("Number of rows in X and Y differ")
#     }

#     if (nrow(Y_fin) != nsamples) {
#         warning("Number of observations != number of samples")
#     }

#     if (nrow(Y2) != nrow(Xmat)) {
#         warning("The full datasets are mismatched")
#     }

#     # Create train and test datasets
#     ntrain <- round(0.7 * nsamples, 0)
#     ntest <- nsamples - ntrain
#     ind_train <- sample(seq_len(nsamples), ntrain)
#     Y_train <- Y_fin[ind_train, ]
#     Yc_train <- Yc_fin[ind_train, ]
#     Y_test <- Y_fin[-ind_train, ]
#     Yc_test <- Yc_fin[-ind_train, ]
#     X_train <- X_fin[ind_train, ]
#     X_test <- X_fin[-ind_train, ]

#     if (binary) {
#         Yr <- Y2
#         Yr_train <- Y_train
#         Yr_test <- Y_test
#     } else {
#         Yr <- Y
#         Yr_train <- Yc_train
#         Yr_test <- Yc_test
#     }

#     if (print) {
#         printDetails(Y_train[, sp], Y_test[, sp], digits = 4L, width = 72L)
#     }

#     e <- list(
#         "Ytrue" = Yr, "X" = Xmat, "Ytrain" = Yr_train, "Xtrain" = X_train,
#         "true_betas" = beta, "train_index" = ind_train,
#         "Ytest" = Yr_test, "Xtest" = X_test
#     )
#     return(e)
# }


simulateAbundanceData_new <- function(grid,
                                      nsamples = 100,
                                      binary = TRUE,
                                      print = TRUE,
                                      tree,
                                      ntraits,
                                      sigma2,
                                      alpha,
                                      epath,
                                      vars = c("CHELSA_bio_1.tif", "CHELSA_bio_13.tif")) { # nolint

    # Get environmental data
    env <- wrangleEnv(vars1 = vars)

    # Blank names
    sp <- tree$tip.label

    #### Simulate responses ####

    ns <- length(sp)

    # +1 for intercept - intercept removed
    # mu <- matrix(0, nrow = (nc), ncol = ns)

    # Expected value
    # mu[1, ] <- 0.2 # Intercept
    # mu[1, ] <- tip_traits$`trait-1`
    # mu[2, ] <- tip_traits$`trait-2`

    # Note vcv assumes BM covariance structure
    # For OU a correlation structure needs to be specified
    # L <- ape::vcv(tree)
    # beta <- apply(mu, 1, function(x) {
    #   MASS::mvrnorm(n = 1, mu = x, Sigma = diag(1, ns, ns))
    # # })

    beta <- simOU_new(
        tree,
        ntraits, sigma2, alpha
    )
    # beta <- cbind(rep(1, ns), rep(1, ns))

    beta <- beta + 3

    # Systematic sampling grid first
    Xmat <- makeXmat(env)
    # Xmat <- Xmat[, -c(which(colnames(Xmat) == "Effort"))]
    M <- Xmat %*% t(beta)
    pool2 <- nrow(Xmat)

    # Independently simulating species observations
    covar <- diag(ns)
    diag(covar) <- 2
    Y <- M + MASS::mvrnorm(n = nrow(M), mu = rep(0, ns), Sigma = covar)

    # Fineaggling
    Y2 <- matrix(as.numeric(Y > 0), ncol = ncol(Y))
    Y2 <- data.frame(Y2)
    Y <- data.frame(Y)
    colnames(Y2) <- colnames(Y) <- sp
    Xmat <- data.frame(Xmat)
    Y2$rwid <- Xmat$rwid <- Y$rwid <- paste0("rw_", seq_len(pool2))

    # Create temporary frame so Y can remain the "groud truth"
    Y_temp <- Y2
    Y_temp$tot <- rowSums(Y_temp[, sp])

    k <- sample(seq_len(nrow(Y_temp)), size = nsamples)
    Y_fin <- Y_temp[k, ]


    # More checks
    X_fin <- Xmat[match(Y_fin$rwid, Xmat$rwid), ]
    Yc_fin <- Y[match(Y_fin$rwid, Y$rwid), ]

    # Another check
    if (any(X_fin$rwid != Y_fin$rwid)) {
        warning("Misalignment in X and Y dataframes")
    }

    if (any(nrow(X_fin) != nrow(Y_fin))) {
        warning("Number of rows in X and Y differ")
    }

    if (nrow(Y_fin) != nsamples) {
        warning("Number of observations != number of samples")
    }

    if (nrow(Y2) != nrow(Xmat)) {
        warning("The full datasets are mismatched")
    }

    # Create train and test datasets
    ntrain <- round(0.7 * nsamples, 0)
    ntest <- nsamples - ntrain
    ind_train <- sample(seq_len(nsamples), ntrain)
    Y_train <- Y_fin[ind_train, ]
    Yc_train <- Yc_fin[ind_train, ]
    Y_test <- Y_fin[-ind_train, ]
    Yc_test <- Yc_fin[-ind_train, ]
    X_train <- X_fin[ind_train, ]
    X_test <- X_fin[-ind_train, ]

    if (binary) {
        Yr <- Y2
        Yr_train <- Y_train
        Yr_test <- Y_test
    } else {
        Yr <- Y
        Yr_train <- Yc_train
        Yr_test <- Yc_test
    }

    if (print) {
        printDetails(Y_train[, sp], Y_test[, sp], digits = 4L, width = 72L)
    }

    e <- list(
        "Ytrue" = Yr, "X" = Xmat, "Ytrain" = Yr_train, "Xtrain" = X_train,
        "true_betas" = beta, "train_index" = ind_train,
        "Ytest" = Yr_test, "Xtest" = X_test
    )
    return(e)
}


# Functions ---
scale_safe <- function(df) {
    numeric_cols <- apply(df, 2, is.numeric) # Identify numeric columns
    rws_complete <- which(complete.cases(df))
    # Apply scaling
    df[rws_complete, numeric_cols] <- apply(df[rws_complete, numeric_cols], 2, function(col) {
        if (length(unique(col)) == 1) {
            return(rep(0, length(col))) # Replace constant columns with 0s
        } else {
            return(scale(col)) # Standard scaling
        }
    })

    return(df)
}

process_x_true_test <- function(x_true_test, cood_test = NULL) {
    x_true_test <- as.data.frame(store_test$x[[fsp1]])
    x_true_test[is.infinite(x_true_test[, "distance"]), "distance"] <- 0
    # Step 1: Add column names
    colnames(x_true_test) <- c(
        "meanTemp", "tempSeasonality", "precipQuart",
        "precipSeasonality", "cloudCover", "Annual_EVI", "TRI", "elevation",
        "distance", "duration", "num_observers"
    )

    # Step 2: Add new predictors
    x_true_test$meanTemp2 <- (x_true_test$meanTemp)^2
    x_true_test$precipQuart2 <- (x_true_test$precipQuart)^2

    # Step 3: Scale the data
    x_true_test <- scale_safe(x_true_test)

    # Step 4: Convert to data frame and add Intercept
    x_true_test <- as.data.frame(x_true_test)
    x_true_test$Intercept <- rep(1, nrow(x_true_test))

    # Step 5: Convert back to matrix
    x_true_test <- as.matrix(x_true_test)

    # Step 6: Reorder columns
    order <- c(
        "Intercept", "meanTemp", "meanTemp2", "tempSeasonality",
        "precipQuart", "precipQuart2", "precipSeasonality", "elevation",
        "Annual_EVI", "TRI", "cloudCover", "duration", "distance",
        "num_observers"
    )
    x_true_test <- x_true_test[, order]

    # Step 7: Remove rows with missing values
    rm_ind <- which(complete.cases(x_true_test) == FALSE)
    if (length(rm_ind) > 1) {
        x_true_test <- x_true_test[-rm_ind, ]
        if (!is.null(cood_test)) {
            cood_test <- cood_test[-rm_ind, ]
        }
    }

    # Return the processed data
    if (is.null(cood_test)) {
        return(x_true_test)
    } else {
        return(list(x_true_test = x_true_test, cood_test = cood_test))
    }
}

calculate_tss <- function(predicted, observed) {
    # Ensure binary classification (0 and 1)
    if (!all(predicted %in% c(0, 1)) | !all(observed %in% c(0, 1))) {
        stop("Predicted and observed values must be binary (0 or 1).")
    }

    # Check if there are at least two unique values in both predicted and observed
    unique_pred <- unique(predicted)
    unique_obs <- unique(observed)

    if (length(unique_pred) == 1 || length(unique_obs) == 1) {
        return(NA) # Not enough class diversity to compute TSS
    }

    # Create a confusion matrix ensuring all possible classes (0,1) exist
    confusion_matrix <- table(
        factor(predicted, levels = c(0, 1)),
        factor(observed, levels = c(0, 1))
    )

    # Extract values safely from confusion matrix
    tp <- confusion_matrix[2, 2] # True Positives
    fn <- confusion_matrix[2, 1] # False Negatives
    fp <- confusion_matrix[1, 2] # False Positives
    tn <- confusion_matrix[1, 1] # True Negatives

    # Avoid division by zero
    sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else NA
    specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA

    # Calculate TSS
    tss <- sensitivity + specificity - 1

    return(tss)
}

gimmeRatios_pp <- function(YFSP, YCOOD, REF_RAST, NAMES = c("x", "y")) {
    if (length(YFSP) != nrow(YCOOD)) {
        error("Something is wrong - Y and COOD are not equal dimensions")
    }
    ind1 <- which(YFSP == 1)
    ind0 <- which(YFSP == 0)
    c0 <- YCOOD[ind0, ]
    c0_vect <- terra::vect(c0[, NAMES], geom = NAMES)
    er0 <- terra::rasterize(c0_vect, REF_RAST, fun = sum)
    er0[is.na(er0)] <- 0
    c1 <- YCOOD[ind1, ]
    c1_vect <- terra::vect(c1[, NAMES], geom = NAMES)
    er1 <- terra::rasterize(c1_vect, REF_RAST, fun = sum)
    er1[is.na(er1)] <- 0
    # err <- (er1 / er0)
    # err[is.infinite(err)] <- 1
    # err[err == 0] <- NA
    p <- terra::as.points(er1)
    data <- cbind(terra::crds(p), terra::values(p))
    return(data)
}

rasterizeCood_pp <- function(YC_CROP, SPS, ER) {
    data <- list()
    for (j in seq_len(length(SPS))) {
        print(SPS[j])
        yfsp <- YC_CROP[, SPS[j]]
        # Get indices for presences/absences
        ind0 <- which(yfsp == 0)
        ind1 <- which(yfsp == 1)
        if (length(ind1) == 0) {
            print("No data")
            data[[j]] <- NA
        } else {
            # Returns the ratio of presences to absences at every pixel
            err_df <- gimmeRatios_pp(yfsp, YC_CROP[, c("x", "y")], ER)
            err_df$sp <- SPS[j]
            data[[j]] <- err_df
        }
    }
    return(data)
}

# Create indices for training and testing sites
# For the Poisson Point Process model
create_indices <- function(STORE, SPLIT) {
    # Initialize a list to store species-specific data
    species_data <- list()
    y_full <- STORE$y
    J <- length(unique(y_full$species))

    # Loop through each species and subset the data
    for (sp in seq_len(J)) {
        species_data[[sp]] <- y_full[y_full$species == sp, ]
    }

    # Initialize lists to store training and testing site indices
    sites_in_training <- list()
    sites_in_testing <- list()

    # Loop through each species-specific data
    for (i in seq_len(length(species_data))) {
        sp_data <- species_data[[i]]
        N <- nrow(sp_data)
        N_train <- floor(N * SPLIT)

        # If the species has less than 3 records, use all for training and none for testing
        if (N < 3) {
            sites_in_training[[i]] <- sp_data$site[1:N]
            sites_in_testing[[i]] <- NA
        } else {
            # Randomly sample indices for training
            idx <- sample(1:N, N_train)
            sites_in_training[[i]] <- sp_data$site[idx]
            sites_in_testing[[i]] <- sp_data$site[-idx]
        }
    }

    # Combine all training and testing site indices
    training_sites <- do.call(c, sites_in_training)
    testing_sites <- do.call(c, sites_in_testing)

    return(list(training_sites = training_sites, testing_sites = testing_sites))
}


# LGCP conditional prediction
# Define a function for conditional prediction
conditional_prediction_LGCP <- function(stan_fit, D_phylo, observed_index, new_index, stan_data) {
    # Extract posterior samples of B from Stan fit
    posterior_samples <- rstan::extract(stan_fit, pars = "B")$B # Shape: (num_samples, J, K)

    # Extract GP hyperparameters from posterior mean estimates
    alpha <- mean(rstan::extract(stan_fit, pars = "alpha")$alpha) # GP amplitude
    rho <- mean(rstan::extract(stan_fit, pars = "rho")$rho) # GP length scale

    # Compute the GP covariance matrix
    K_phylo <- alpha^2 * exp(-D_phylo^2 / (2 * rho^2))

    # Partition the covariance matrix using numeric indices
    K_obs_obs <- K_phylo[observed_index, observed_index] # Covariance among observed species
    K_new_obs <- K_phylo[new_index, observed_index, drop = FALSE] # Covariance between new and observed species
    K_new_new <- K_phylo[new_index, new_index, drop = FALSE] # Variance for the new species

    # Compute posterior mean for observed B values
    B_obs_mean <- apply(posterior_samples[, observed_index, ], c(2, 3), mean) # Posterior mean for B_obs

    # Compute conditional mean (convert to vector)
    mu_new <- as.numeric(K_new_obs %*% solve(K_obs_obs, B_obs_mean)) # Shape: (K,)

    # Compute conditional covariance
    Sigma_new <- K_new_new - K_new_obs %*% solve(K_obs_obs, t(K_new_obs))

    # Ensure Sigma_new is a diagonal covariance matrix for sampling
    if (length(new_index) == 1) {
        sigma_new <- sqrt(as.numeric(Sigma_new)) # Convert scalar variance to SD
    } else {
        sigma_new <- sqrt(diag(Sigma_new)) # Extract standard deviations for each predictor
    }

    # Initialize matrix to store sampled B values
    num_samples <- dim(posterior_samples)[1] # Number of MCMC samples
    B_new_samples <- matrix(NA, nrow = num_samples, ncol = K)

    # Loop over each environmental variable and sample from univariate normal
    for (k in seq_len(K)) {
        B_new_samples[, k] <- rnorm(num_samples, mean = mu_new[k], sd = sigma_new)
    }

    # Convert to dataframe
    B_new_df <- data.frame(B_new_samples)

    colnames(B_new_df) <- colnames(stan_data$X)

    # Return the inferred B values as a dataframe
    return(B_new_df)
}


# conditional_prediction_LGCP <- function(stan_fit, D_phylo, observed_index, new_index, stan_data,
#                                         st_devs = NULL, jitter = 1e-6, use_draws = TRUE) {
#     B_draws <- rstan::extract(stan_fit, pars = "B")$B # [S, J, K]
#     alpha_draws <- rstan::extract(stan_fit, pars = "alpha")$alpha
#     rho_draws <- rstan::extract(stan_fit, pars = "rho")$rho

#     S <- dim(B_draws)[1]
#     J <- dim(B_draws)[2]
#     K <- dim(B_draws)[3]

#     if (!use_draws) { # optional plug-in mode (faster but less correct)
#         B_draws <- array(apply(B_draws, c(2, 3), mean), dim = c(1, J, K))
#         alpha_draws <- mean(alpha_draws)
#         rho_draws <- mean(rho_draws)
#         S <- 1
#     }

#     ker <- function(a, r, D) a^2 * exp(-(D^2) / (2 * r^2))

#     out <- matrix(NA_real_, nrow = S, ncol = K)

#     for (s in 1:S) {
#         Kfull <- ker(alpha_draws[s], rho_draws[s], D_phylo)

#         K_RR <- Kfull[observed_index, observed_index, drop = FALSE]
#         K_DR <- Kfull[new_index, observed_index, drop = FALSE]
#         K_DD <- Kfull[new_index, new_index, drop = FALSE] # scalar if single new

#         if (!is.null(st_devs)) K_RR <- K_RR + diag(st_devs[observed_index]^2)
#         L <- chol(K_RR + diag(jitter, nrow(K_RR)))
#         solve_KRR <- function(v) backsolve(L, forwardsolve(t(L), v))

#         for (k in 1:K) {
#             B_R_k <- B_draws[s, observed_index, k]
#             mu_k <- as.numeric(K_DR %*% solve_KRR(B_R_k))
#             var_k <- as.numeric(K_DD - K_DR %*% solve_KRR(t(K_DR)))
#             sd_k <- sqrt(max(var_k, jitter))
#             out[s, k] <- rnorm(1, mu_k, sd_k) # or store mu_k if you don't want sampling
#         }
#     }

#     df <- as.data.frame(out)
#     df
# }


# Simple, stable per-draw conditional predictor (species-level GP over betas)
# Returns: data.frame with S rows (posterior draws) and K columns (predictors)
# Notes:
# - Assumes new_index is a single focal species.
# - Variance is species-level (same sd for all predictors), which matches your model.
# - Optional nugget on K_RR if you have species-specific st_devs draws.

cp_LGCP_simple <- function(
    stan_fit,
    D_phylo, # J x J distance matrix among species
    observed_index, # integer vector of observed-species indices (R)
    new_index, # single integer index of focal species (D)
    stan_data, # to pick K and column names
    use_draws = TRUE, # TRUE: propagate uncertainty; FALSE: plug-in means (S=1)
    st_devs_draws = NULL, # optional S x J matrix of species SDs; if provided, adds nugget
    jitter = 1e-6 # numeric jitter for numerical stability
    ) {
    stopifnot(length(new_index) == 1)

    # --- Extract posterior draws ---
    B_draws <- rstan::extract(stan_fit, pars = "B")$B # [S, J, K]
    alpha_draws <- rstan::extract(stan_fit, pars = "alpha")$alpha
    rho_draws <- rstan::extract(stan_fit, pars = "rho")$rho

    S <- dim(B_draws)[1]
    J <- dim(B_draws)[2]
    K <- dim(B_draws)[3]

    # Optional plug-in (fast, less correct): collapse to means
    if (!use_draws) {
        B_draws <- array(apply(B_draws, c(2, 3), mean), dim = c(1, J, K))
        alpha_draws <- mean(alpha_draws)
        rho_draws <- mean(rho_draws)
        if (!is.null(st_devs_draws)) {
            # collapse st_devs if provided as draws
            st_devs_draws <- matrix(colMeans(st_devs_draws), nrow = 1) # [1, J]
        }
        S <- 1
    }

    # Kernel function (squared exponential on phylo distances)
    ker <- function(a, r, D) a^2 * exp(-(D^2) / (2 * r^2))

    out <- matrix(NA_real_, nrow = S, ncol = K)

    for (s in seq_len(S)) {
        # Build kernel for this draw
        Kfull <- ker(alpha_draws[s], rho_draws[s], D_phylo)

        # Partition by species
        K_RR <- Kfull[observed_index, observed_index, drop = FALSE]
        K_DR <- Kfull[new_index, observed_index, drop = FALSE] # 1 x |R|
        K_DD <- Kfull[new_index, new_index, drop = FALSE] # scalar

        # Optional nugget on observed block (species-specific residual sd)
        if (!is.null(st_devs_draws)) {
            stopifnot(ncol(st_devs_draws) == J)
            K_RR <- K_RR + diag(st_devs_draws[s, observed_index]^2, nrow = length(observed_index))
        }

        # Stable solve via Cholesky + jitter
        L_RR <- chol(K_RR + diag(jitter, nrow(K_RR)))
        solve_KRR <- function(B) backsolve(L_RR, forwardsolve(t(L_RR), B))

        # Collect observed betas for ALL predictors at once: (|R| x K)
        B_R_mat <- do.call(cbind, lapply(seq_len(K), function(k) B_draws[s, observed_index, k]))

        # Conditional mean for ALL predictors: (1 x |R|) %*% (|R| x K) = (1 x K)
        mu_vec <- as.vector(K_DR %*% solve_KRR(B_R_mat))

        # Conditional variance is species-level (same for each predictor)
        var_D <- as.numeric(K_DD - K_DR %*% solve_KRR(t(K_DR)))
        sd_D <- sqrt(max(var_D, jitter))

        # One posterior-predictive draw per predictor for this draw s
        out[s, ] <- rnorm(K, mean = mu_vec, sd = sd_D)
    }

    df <- as.data.frame(out)
    colnames(df) <- colnames(stan_data$X)
    df
}
