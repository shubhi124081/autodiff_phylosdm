#include <TMB.hpp>

template <class Type>
Type objective_function<Type>::operator()()
{
  using namespace density;

  // ------------------ DATA ------------------
  DATA_INTEGER(N);
  DATA_INTEGER(J);
  DATA_INTEGER(K);
  DATA_IVECTOR(species); // 1-based from R
  DATA_MATRIX(X);        // N x K
  DATA_VECTOR(y);        // N (0/1 here but Poisson)
  DATA_MATRIX(D_phylo);  // J x J
  DATA_VECTOR(offset);   // N

  // ---------------- PARAMETERS ---------------
  PARAMETER_MATRIX(B); // J x K

  // positive params via log-transform
  PARAMETER(log_alpha);
  PARAMETER(log_rho);
  PARAMETER(log_sigma_f);

  PARAMETER_VECTOR(f); // N

  // transforms
  Type alpha = exp(log_alpha);
  Type rho = exp(log_rho);
  Type sigma_f = exp(log_sigma_f);

  Type nll = 0.0;

  // Priors on original (positive) scale + Jacobians
  nll -= dnorm(alpha, Type(0), Type(1), true);
  nll -= log(alpha);
  nll -= dnorm(rho, Type(0), Type(1), true);
  nll -= log(rho);
  nll -= dnorm(sigma_f, Type(0), Type(1), true);
  nll -= log(sigma_f);

  // ---------- phylogenetic covariance (RBF) ----------
  matrix<Type> Kphy(J, J);
  for (int i = 0; i < J; ++i)
  {
    for (int j = 0; j < J; ++j)
    {
      Type d = D_phylo(i, j);
      Kphy(i, j) = alpha * alpha * exp(-(d * d) / (Type(2.0) * rho * rho));
    }
  }
  for (int j = 0; j < J; ++j)
    Kphy(j, j) += Type(1e-6); // jitter

  // MVN prior on columns of B
  MVNORM_t<Type> mvn(Kphy);
  for (int k = 0; k < K; ++k)
  {
    vector<Type> Bk = B.col(k);
    nll += mvn(Bk);
  }

  // Spatial RE prior: f ~ N(0, sigma_f^2 I)
  for (int n = 0; n < N; ++n)
  {
    nll -= dnorm(f(n), Type(0), sigma_f, true);
  }

  // Likelihood: Poisson with log link
  for (int n = 0; n < N; ++n)
  {
    int s = species(n) - 1; // 0-based

    // ---- FIX: robust dot product without Eigen AD gymnastics ----
    Type xb = 0.0;
    for (int k = 0; k < K; ++k)
    {
      xb += X(n, k) * B(s, k);
    }
    // -------------------------------------------------------------

    Type log_lambda = xb + f(n) + offset(n);
    nll -= dpois(y(n), exp(log_lambda), true);
  }

  ADREPORT(alpha);
  ADREPORT(rho);
  ADREPORT(sigma_f);

  return nll;
}
