data {
  int<lower=0> n;                 // Number of observations
  vector[n] V;                    // Confluency variable (response variable)
  int<lower=1> K;                 // Number of groups (treatments)
  int<lower=1, upper=K> group[n]; // Group assignment for each observation
}

parameters {
  vector<lower=1>[K] a;           // Shape parameter 'a' for each group (α > 1 to ensure concave, unimodal Beta distribution)
  vector<lower=1>[K] b;           // Shape parameter 'b' for each group (β > 1 to ensure concave, unimodal Beta distribution)
  vector<lower=0>[K] zeta;
  vector<lower=0>[K] varrho;
  vector<lower=0>[K] delta;
  vector<lower=0>[K] xi;

}

model {
  // Priors for hyperparameters (zeta, varrho, delta, xi)

  for (k in 1:K) {

    zeta[k] ~ gamma(0.5, 1);
   varrho[k] ~ gamma(2, 2);
     delta[k] ~ gamma(5, 1);
   xi[k] ~ gamma(1, 1);
  a[k] ~ gamma(zeta[k], varrho[k]);       // Prior
  b[k] ~ gamma(delta[k], xi[k]);      // Prior
  }

  // Likelihood: Beta distribution with concave, unimodal shape for confluency variable V
  for (j in 1:n) {
    target += beta_lpdf(V[j] | a[group[j]], b[group[j]]);  // Log-likelihood for each observation
  }
}

generated quantities {
  vector[n] ypred;    // Predicted values from posterior predictive distribution
  vector[n] log_lik;  // Log-likelihood for each observation

  // Posterior predictive sampling and log-likelihood calculation
  for (j in 1:n) {
    ypred[j] = beta_rng(a[group[j]], b[group[j]]);         // Generate predictive samples
    log_lik[j] = beta_lpdf(V[j] | a[group[j]], b[group[j]]); // Compute log-likelihood
  }
}

