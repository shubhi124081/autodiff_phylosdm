library(TMB)

# --- toy data just to run ---
set.seed(1)
N <- 1000
J <- 48
K <- 2
species <- sample(1:J, N, replace = TRUE)

# True parms for sim (simple)
B_true <- matrix(1, J, K) # species-by-covariate
sigma_f_true <- 0.5
f_true <- rnorm(N, 0, sigma_f_true) # iid site RE (matches your model)
offset <- rep(0, N)
X <- matrix(rnorm(N * K), N, K)
D_phylo <- as.matrix(dist(1:J))

# Per-observation linear predictor: X[n,]*B[species[n],]
eta <- rowSums(X * B_true[species, , drop = FALSE]) + f_true + offset
lambda <- exp(eta) # Poisson rate must be positive
y <- rpois(N, lambda)

# data <- list(
#     N = N, J = J, K = K, species = species, X = X,
#     y = y, D_phylo = D_phylo, offset = offset
# )

# compile("lgcp_background.cpp")
# dyn.load(dynlib("lgcp_background"))


# par <- list(
#     B = matrix(0, J, K),
#     log_alpha = log(0.5),
#     log_rho = log(1.0),
#     log_sigma_f = log(0.5),
#     f = rep(0, N)
# )

# S1 <- Sys.time()
# obj <- MakeADFun(data, par, DLL = "lgcp_background", random = c("f", "B"))
# opt <- nlminb(obj$par, obj$fn, obj$gr)
# rep <- sdreport(obj)
# summary(rep, "report") # alpha, rho, sigma_f
# S2 <- Sys.time()
# mins <- as.numeric(difftime(S2, S1, units = "mins"))


# # ---- Get B estimates (since B is 'random') ----
# rnd <- summary(rep, "random")[, 1] # first column = means
# # Ordering of 'random' matches the concatenation of parameters you listed: vec(B), then f
# B_hat <- matrix(rnd[seq_len(J * K)], nrow = J, ncol = K, byrow = FALSE)

# cat("\nColumn means of B_hat:\n")
# print(colMeans(B_hat))

# Run stan model
library(rstan)
source("~/TMB_LGCP/stan.R")
code <- as.character(LGCP_background)
S1 <- Sys.time()
mod <- rstan::stan(
    model_code = code,
    data = list(
        N = N,
        J = J,
        K = K,
        E = 1,
        X = X,
        y = y,
        D_phylo = D_phylo,
        offset = offset,
        species = species
    ),
    chains = 1,
    iter = 8000,
    warmup = 3500,
    control = list(
        adapt_delta = 0.9,
        stepsize = 1,
        stepsize_jitter = 0
    )
)
S2 <- Sys.time()
mins <- as.numeric(difftime(S2, S1, units = "mins"))
cat("Stan time (mins):", mins, "\n")


# fit <- as.matrix(mod)
# colMeans(fit[, grep("B", colnames(fit))])
# mean(fit[, grep("alpha", colnames(fit))])
# mean(fit[, grep("rho", colnames(fit))])

save(mod, file = "~/Downloads/stan_fit_J48.RData")
