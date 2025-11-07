# Set-up --
# This script is currently interactive!
HPC <- Sys.getenv("HPC")
if (HPC == "FALSE") {
    root <- "~/phylo-sdms2"
} else {
    root <- "~/phylo-sdms2"
}

data_directory <- file.path(root, "data")
scripts_directory <- file.path(root, "scripts")
res_directory <- file.path(root, "res")
source(file.path(scripts_directory, "000-phyloGenie_functions.R"))

# Define dataset and model arguments
EXP_ROOT <- "happy" # folder in analysis
EXP_ID <- "temp" # exp_root + exp_id = experiment_name
CLUSTER <- "Phaethornis2" # from spList, exp_root + exp_id + cluster = shortname
FSP <- "ALL" # ALL if all species otherwise species name "GENUS_SPECIES"
DEF_LEV <- "t2008" # e2024 is FULL
REPNO <- "1" # Rep number, not number of reps
MODEL_NAME <- "LGCP_background" # must be from 012-stan_models.R
# Alternatively, could call config file

# Load model
shortname <- paste0(EXP_ROOT, "_", EXP_ID, "_", CLUSTER)
data_name <- paste0(shortname, "_", FSP, "_", DEF_LEV, "_", REPNO)
shortname_model <- paste0(data_name, "_", MODEL_NAME)
contents <- load(file.path(res_directory, shortname, paste0(shortname_model, ".Rdata"))) # loads object all

# 1. Print parameter summary + Rhat
config <- all$config
fit_lgcp <- all$fit
print(fit_lgcp, probs = c(0.025, 0.5, 0.975), digits = 2)

# 2. Check convergence (Rhat < 1.01 ideally)
summary_fit <- rstan::summary(fit_lgcp)$summary
rhat_vals <- summary_fit[, "Rhat"]
cat("Parameters with Rhat > 1.01:\n")
print(names(rhat_vals[rhat_vals > 1.01]))

# 3. Traceplots for key parameters
# Before plots, open a plot window
httpgd::hgd()
httpgd::hgd_browse()
# Choose a few parameters (e.g., beta coefficients and spatial effects)
rstan::traceplot(fit_lgcp, pars = c("alpha", "rho", "sigma_f"), inc_warmup = FALSE)

# For species-environment coefficients:
# Traceplot for B[1,1], B[2,1], ..., as a sample
rstan::traceplot(fit_lgcp, pars = c("B[1,1]", "B[2,1]", "B[3,1]"))

# 4. Posterior Predictive Check (simple)
# Posterior predictive samples of y
# Extract posterior samples
B_post <- rstan::extract(fit_lgcp, pars = "B")$B # iterations × J × K
f_post <- rstan::extract(fit_lgcp, pars = "f")$f # iterations × N_obs

# Data inputs
# Load data used for fitting
contents <- load(file.path(data_directory, shortname, paste0(data_name, ".Rdata")))
standata <- everything$data # data used for fitting
X <- standata$X # matrix[N_obs, K]
species <- standata$species # vector[N_obs]
offset <- standata$offset # vector[N_obs]
y_obs <- standata$y # observed counts or presences

# Set up matrix for posterior predictive samples
n_iter <- dim(f_post)[1]
n_obs <- length(y_obs)
lambda_post <- matrix(NA, nrow = n_iter, ncol = n_obs)

# Loop to compute posterior predicted λ for each observation
for (i in 1:n_iter) {
    for (n in 1:n_obs) {
        s <- species[n]
        xb <- sum(X[n, ] * B_post[i, s, ])
        eta <- xb + f_post[i, n] + log(abs(offset))[n]
        lambda_post[i, n] <- rpois(1, exp(eta))
    }
}

# Posterior predictive mean
ppc_mean <- colMeans(lambda_post)

# Plot observed vs predicted
hist(ppc_mean,
    breaks = 30,
    main = "Predicted intensity (λ) at observed presences",
    xlab = "Predicted λ", col = "skyblue"
)

species_ids <- standata$species
boxplot(split(ppc_mean, species_ids),
    main = "Predicted λ by species",
    xlab = "Species", ylab = "Predicted λ", col = "lightgreen"
)

# Add uncertainty bands
pred_lower <- apply(lambda_post, 2, quantile, probs = 0.025)
pred_upper <- apply(lambda_post, 2, quantile, probs = 0.975)
arrows(x0 = y_obs, y0 = pred_lower, y1 = pred_upper, angle = 90, code = 3, length = 0.02, col = "gray60")

# Write out plots
# MCMC evaluation script
print("Running 030-MCMC_eval.R ....")

# Load result file ----
files <- list.files(file.path(res_directory, shortname))

p <- list()
dict <- data.frame(
    "variable" = c("rhats", "ess_bulks", "ess_tails"),
    "desired_min" = c(1.05, 100, 100)
)

for (i in seq_len(length(files))) {
    contents <- load(file.path(root, "res", exp_name, files[i]))
    fit <- as.matrix(all$fit)
    rhats <- apply(fit, 2, rstan::Rhat)
    ess_bulks <- apply(fit, 2, rstan::ess_bulk)
    ess_tails <- apply(fit, 2, rstan::ess_tail)

    stat <- data.frame(
        "rhats" = rhats,
        "ess_bulks" = ess_bulks,
        "ess_tails" = ess_tails
    )

    stat <- reshape2::melt(stat)

    p[[i]] <- ggplot(stat) +
        geom_histogram(aes(x = value)) +
        geom_vline(data = dict, aes(xintercept = desired_min), col = "red") +
        facet_wrap("~variable", scales = "free") +
        theme_bw() +
        ggtitle(files[i])
}

p_arr <- gridExtra::marrangeGrob(p, nrow = 1, ncol = 1)

# Save arrangement -----
ggsave(file.path(
    root, "analysis", exp_root, folder, "figures",
    paste0(exp_name, "_mcmc_summary.pdf")
), p_arr)
