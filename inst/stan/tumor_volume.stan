data {
  int<lower=0> n;                   // Number of observations
  vector[n] V;                      // Tumor volume variable (response variable)
  int<lower=1> K;                   // Number of groups (treatments)
  int<lower=1, upper=K> group[n];   // Group assignment for each observation
}

parameters {
  vector[K] mu;                     // Mean for each group
  vector<lower=0>[K] sigma;         // Standard deviation for each group
  vector<lower=0>[K] nu;                 // Degrees of freedom for the Student-t distribution
  vector<lower=0>[K] sigma_raw;          // Scale parameter for the group means (mu)
  vector[K] mu_mu;              // Hyperparameter mean for group means (mu)
  vector<lower=0>[K] rate;               // Rate parameter for the Exponential distribution on sigma
}

model {


  // Priors for each group's mean and standard deviation
  for (k in 1:K) {
      // Hyperprior for the global mean (mu_mu)
  mu_mu[k] ~ normal(0, 5);

  // Prior for the degrees of freedom of the Student-t distribution
  nu[k] ~ exponential(1);

  // Prior for the scale parameter of the group means (mu) with a Gamma distribution
  sigma_raw[k] ~ gamma(4, 2);

  // Prior for the rate parameter for the Exponential distribution on sigma
  rate[k] ~ gamma(4, 2);
    mu[k] ~ student_t(nu[k], mu_mu[k], sigma_raw[k]);  // Student-t prior for group means
    sigma[k] ~ exponential(rate[k]);             // Exponential prior for group standard deviations
  }

  // Likelihood function for each observation
  for (i in 1:n) {
    target += lognormal_lpdf(V[i] | mu[group[i]], sigma[group[i]]);  // Lognormal likelihood
  }
}

generated quantities {
  vector[n] ypred;       // Posterior predictive samples for each observation
  vector[n] log_lik;     // Log-likelihood for each observation

  // Loop to generate predictions and calculate log-likelihoods
  for (i in 1:n) {
    ypred[i] = lognormal_rng(mu[group[i]], sigma[group[i]]);  // Generate posterior predictive values
    log_lik[i] = lognormal_lpdf(V[i] | mu[group[i]], sigma[group[i]]);  // Calculate log-likelihood
  }
}
