library(RTMB)

set.seed(41)

# Dimensions
N <- 200 # Obsersvations
J <- 1 # No. of species
K <- 6 # No. of covariates

# We're just going to simulate a single species poisson model
species <- sample.int(J, N, replace = TRUE)
