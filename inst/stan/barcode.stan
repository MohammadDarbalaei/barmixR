// -----------------------------------------------------
// define dirichlet_multinomial for likelihood estimation
// -----------------------------------------------------
functions {

  real dirichlet_multinomial_lpmf(int[] y, vector alpha) {
    real alpha_plus = sum(alpha);
    return lgamma(alpha_plus) + sum(lgamma(alpha + to_vector(y)))
                - lgamma(alpha_plus+sum(y)) - sum(lgamma(alpha));
  }

  int[] dirichlet_multinomial_rng(vector alpha, int N) {
    return multinomial_rng(dirichlet_rng(alpha), N);
  }

}
// -----------------------------------------------------
// -----------------------------------------------------

data {
  int<lower=0> n; // Number of observations (pools or tumor)
  int<lower=0> C; // cell lines * (replicate of each cell line) >> for example 8*3
  int y [n,C]; // observed cell counts
  int<lower=1> K;                   // Number of groups (number of treatments)
  int<lower=1, upper=K> group[n];    // Group assignment for each observation

  // Optional hyperparameters for psi and varphi priors
  real<lower=0> psi_mean; // Mean of psi's lognormal prior
  real<lower=0> psi_sd;   // Standard deviation of psi's lognormal prior
  real<lower=0> varphi_mean; // Mean of varphi's lognormal prior
  real<lower=0> varphi_sd;   // Standard deviation of varphi's lognormal prior
}

transformed data {
  real default_psi_mean = 0;
  real default_psi_sd = 0.5;
  real default_varphi_mean = 0;
  real default_varphi_sd = 0.01;

  real final_psi_mean = is_nan(psi_mean) ? default_psi_mean : psi_mean;
  real final_psi_sd = is_nan(psi_sd) ? default_psi_sd : psi_sd;
  real final_varphi_mean = is_nan(varphi_mean) ? default_varphi_mean : varphi_mean;
  real final_varphi_sd = is_nan(varphi_sd) ? default_varphi_sd : varphi_sd;
}

parameters {
  matrix[K,C] beta;   // treatment scores
  vector<lower = 0>[C] S;  // precision
  vector<lower = 0>[C] psi; // mean of S (precision)
  vector<lower = 0>[C] varphi; // standard deviation of S (precision)
  vector[C] w; // mean of mu (hyperparameter)
  matrix[K,C] mu; // matrix of observation means
  real<lower=0> tau; // between categories
  matrix<lower=0>[K, C] sigma; // within categories (updated to 2D)
  vector<lower=0>[K] nu; // for Student-t noise
}

transformed parameters {
  simplex[C] theta[n];
  vector[C] x_beta[n];

  for (j in 1:n) {
    x_beta[j] = to_vector(beta[group[j],]);
    theta[j] = softmax(x_beta[j]);
  }
}

model {
  // Priors
  psi ~ lognormal(final_psi_mean, final_psi_sd);
  varphi ~ lognormal(final_varphi_mean, final_varphi_sd);
  S ~ lognormal(psi, varphi);
  w ~ normal(0, 1);
  tau ~ gumbel(0.5, 0.1);
  nu ~ exponential(1);

  for (k in 1:K) {
    for (l in 1:C) {
      mu[k, l] ~ normal(w[l], tau);
      sigma[k, l] ~ gumbel(0.5, 0.1);
    }
    beta[k, ] ~ student_t(nu[k], mu[k, ], sigma[k, ]);
  }

  for (j in 1:n) {
    y[j] ~ dirichlet_multinomial(S .* theta[j]);
  }
}

generated quantities {
  int y_rep[n, C];
  real log_lik[n];

  // Simulations
  for (j in 1:n) {
    y_rep[j] = dirichlet_multinomial_rng(S .* theta[j], sum(y[j]));
    log_lik[j] = dirichlet_multinomial_lpmf(y[j] | S .* theta[j]);
  }
}
