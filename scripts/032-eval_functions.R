predict_intensity_rf <- function(posterior_B, posterior_sigma_f, x_new, y_new, offset_new) {
    # Extract dimensions
    num_samples <- dim(posterior_B)[1]
    N_new <- nrow(x_new)

    log_lambda_samples <- matrix(NA, nrow = num_samples, ncol = N_new)
    pb <- txtProgressBar(min = 0, max = num_samples, style = 3)

    for (i in 1:num_samples) {
        for (n in 1:N_new) {
            species_n <- y_new$species[n] # species index
            f_new <- rnorm(1, mean = 0, sd = posterior_sigma_f[i])
            log_lambda_samples[i, n] <- sum(x_new[n, ] * posterior_B[i, species_n, ]) +
                f_new + offset_new[n]
        }
        setTxtProgressBar(pb, i)
    }
    close(pb)

    # Return mean predicted intensity (lambda)
    lambda_samples <- exp(log_lambda_samples)
    predicted_intensity <- colMeans(lambda_samples)

    return(predicted_intensity)
}

process_offset <- function(x_test, y_test, offset_vars) {
    offset_test <- x_test[, offset_vars, drop = FALSE]
    offset_test <- cbind(offset_test, site = y_test$site)
    offset_matrix <- offset_test
    colnames(offset_matrix) <- c(offset_vars, "site")

    # Multiply all offset_vars together for each row
    offset_multiplied <- apply(offset_matrix[, offset_vars, drop = FALSE], 1, prod)
    offset_multiplied <- cbind(offset_multiplied, offset_matrix[, "site"])
    colnames(offset_multiplied) <- c("offset", "site")

    offset_test <- abs(offset_multiplied[, "offset"])
    offset_test <- log(offset_test)
    x_test <- x_test[, setdiff(colnames(x_test), offset_vars), drop = FALSE]
    return(list(offset_test = offset_test, x_test = x_test))
}
