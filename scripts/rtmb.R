## ============================================================
## RTMB LGCP-style model with phylogenetic RBF prior on B
## - Robust to common sdreport() issues
## - Non-centered parameterization for B (B = L %*% Z)
## - Priors on LOG scale for (alpha, rho, sigma_f)
## ============================================================

## ---------- 0) Packages ----------
pkgs <- c("Matrix", "RTMB")
to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
lapply(pkgs, library, character.only = TRUE)

## ---------- 1) Simulate data ----------
set.seed(42)

## Dimensions (tweak freely)
N <- 200 # observations
J <- 12 # species
K <- 6 # covariates

## True hyperparameters
alpha_true <- 0.7 # magnitude of phylo kernel
rho_true <- 1.2 # length-scale of phylo kernel
sigma_f_true <- 0.5 # SD of random effect f

## Design + species mapping
species <- sample.int(J, N, replace = TRUE)
X <- matrix(rnorm(N * K), N, K)
X <- scale(X) # helps conditioning
offset <- rep(0, N)

## Simple phylogenetic "distance" (plug your real D_phylo here)
D_phylo <- as.matrix(dist(seq_len(J)))

## RBF kernel for truth
Kphy_true <- alpha_true^2 * exp(-(D_phylo^2) / (2 * rho_true^2))
diag(Kphy_true) <- diag(Kphy_true) + 1e-4
L_true <- chol(Kphy_true)

## True B via non-centered param (Z ~ N(0,1))
Z_true <- matrix(rnorm(J * K), J, K)
B_true <- L_true %*% Z_true # J x K

## Random effect f
f_true <- rnorm(N, 0, sigma_f_true)

## Linear predictor and counts
eta <- rowSums(X * B_true[species, , drop = FALSE]) + f_true + offset
mu <- exp(eta)
y <- rpois(N, lambda = mu)

## Pack data list
data <- list(
    N = N, J = J, K = K,
    species = species,
    X = X,
    y = y,
    D_phylo = D_phylo,
    offset = offset
)

## ---------- 2) RTMB model ----------
rtmb_model <- function(parms) {
    with(c(parms), {
        nll <- 0

        ## Transforms
        alpha <- exp(log_alpha)
        rho <- exp(log_rho)
        sigma_f <- exp(log_sigma_f)

        ## Priors on LOG-scale (no Jacobians)
        ## Centered loosely near the true values; widen if you like
        nll <- nll - dnorm(log_alpha, mean = log(0.7), sd = 0.5, log = TRUE)
        nll <- nll - dnorm(log_rho, mean = log(1.2), sd = 0.5, log = TRUE)
        nll <- nll - dnorm(log_sigma_f, mean = log(0.5), sd = 0.5, log = TRUE)

        ## Phylogenetic RBF kernel and its Cholesky
        Kphy <- alpha^2 * exp(-(D_phylo^2) / (2 * rho^2))
        diag(Kphy) <- diag(Kphy) + 1e-4 # jitter helps PD
        L <- chol(Kphy) # upper-triangular

        ## Non-centered prior: B = L %*% Z, with Z ~ N(0, I)
        ## Z is J x K parameter matrix (fixed)
        nll <- nll - sum(dnorm(Z, 0, 1, log = TRUE))
        B <- L %*% Z # J x K (do NOT transpose)

        ## Random effect: f ~ N(0, sigma_f^2 I)
        nll <- nll - sum(dnorm(f, 0, sigma_f, log = TRUE))

        ## Poisson log-likelihood
        xb <- rowSums(X * B[species, , drop = FALSE])
        eta <- xb + f + offset
        nll <- nll - sum(dpois(y, exp(eta), log = TRUE))

        ## Report transformed params for easy access
        ADREPORT(alpha)
        ADREPORT(rho)
        ADREPORT(sigma_f)

        nll
    })
}

## ---------- 3) Initial values ----------
par <- list(
    log_alpha = log(0.5),
    log_rho = log(1.0),
    log_sigma_f = log(0.5),
    Z = matrix(0, J, K), # non-centered coefficients
    f = rep(0, N) # random effects
)

## ---------- 4) Build, optimize ----------
obj <- RTMB::MakeADFun(
    f = rtmb_model,
    parameters = par,
    random = "f", # integrate out f via Laplace
    silent = TRUE
)

opt <- nlminb(
    start = obj$par,
    objective = obj$fn,
    gradient = obj$gr,
    control = list(eval.max = 5000, iter.max = 5000, rel.tol = 1e-10)
)

cat("Convergence code (0=OK):", opt$convergence, "\n")
cat("Final objective:", opt$objective, "\n\n")

## ---------- 5) Standard errors / reports ----------
rep <- try(RTMB::sdreport(obj), silent = TRUE)
if (inherits(rep, "try-error")) {
    message("sdreport() failed; retrying with skip.delta.method = TRUE (no SEs)")
    rep <- RTMB::sdreport(obj, skip.delta.method = TRUE)
}

## Pull reported quantities (alpha, rho, sigma_f)
sum_rep <- try(summary(rep, "report"), silent = TRUE)

## ---------- 6) Reconstruct B_hat, compare to truth ----------
## Get estimated parameters at optimum
pl <- obj$env$parList(opt$par)
alpha_hat <- exp(pl$log_alpha)
rho_hat <- exp(pl$log_rho)
sigmaf_hat <- exp(pl$log_sigma_f)
Z_hat <- matrix(pl$Z, J, K)

Kphy_hat <- alpha_hat^2 * exp(-(D_phylo^2) / (2 * rho_hat^2))
diag(Kphy_hat) <- diag(Kphy_hat) + 1e-4
L_hat <- chol(Kphy_hat)
B_hat <- L_hat %*% Z_hat

## ---------- 7) Print summaries ----------
cat("--- Truth vs Estimates (report) ---\n")
if (!inherits(sum_rep, "try-error")) {
    print(sum_rep)
    est <- setNames(sum_rep[, "Estimate"], rownames(sum_rep))
    se <- sum_rep[, "Std. Error"]
    cat(sprintf(
        "alpha:   true = %.3f; est = %.3f (SE = %s)\n",
        alpha_true, est["alpha"], ifelse(is.finite(se[rownames(sum_rep) == "alpha"]),
            sprintf("%.3f", se[rownames(sum_rep) == "alpha"]), "NA"
        )
    ))
    cat(sprintf(
        "rho:     true = %.3f; est = %.3f (SE = %s)\n",
        rho_true, est["rho"], ifelse(is.finite(se[rownames(sum_rep) == "rho"]),
            sprintf("%.3f", se[rownames(sum_rep) == "rho"]), "NA"
        )
    ))
    cat(sprintf(
        "sigma_f: true = %.3f; est = %.3f (SE = %s)\n\n",
        sigma_f_true, est["sigma_f"], ifelse(is.finite(se[rownames(sum_rep) == "sigma_f"]),
            sprintf("%.3f", se[rownames(sum_rep) == "sigma_f"]), "NA"
        )
    ))
} else {
    cat("(No SEs available; sdreport fell back or failed)\n")
    cat(sprintf("alpha:   true = %.3f; est = %.3f\n", alpha_true, alpha_hat))
    cat(sprintf("rho:     true = %.3f; est = %.3f\n", rho_true, rho_hat))
    cat(sprintf("sigma_f: true = %.3f; est = %.3f\n\n", sigma_f_true, sigmaf_hat))
}

## Correlation between vec(B_true) and vec(B_hat)
corr_B <- cor(as.vector(B_true), as.vector(B_hat))
cat(sprintf("corr(vec(B_true), vec(B_hat)) = %.3f\n", corr_B))

## Optional quick diagnostics
maxgrad <- max(abs(obj$gr(opt$par)))
cat(sprintf("max |grad| at optimum: %.3e\n", maxgrad))
