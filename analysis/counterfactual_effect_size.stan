functions {
  // given the necessity strength n, the sufficiency strength s,
  // and the prior probability p, calculate the causal strength
  // of a variable
  real K_NS(real n, real s, real p) {
    return (1-p)*n + p*s;
  }
}

data {
  int<lower=1> N;                           // number of data points
  int<lower=1> V;                           // number of vignettes (2: reversal/adversarial)
  int<lower=1> F;                           // number of factors being judged (2: PC/DP)
  array[N] int<lower=1, upper=V> vignette;  // vignette per observation (1: reversal, 2: adversarial)
  array[N] int<lower=1, upper=F> factor;    // factor per observation (1: PC, 2: DP)
  vector<lower=0, upper=1>[N] cause;        // causal judgments
  int<lower=0, upper=1> prior_only;         // sample from just the prior (1) or the posterior (0)?
}

parameters {
  real<lower=0, upper=1> p_PC;
  vector<lower=0, upper=1>[V] p_PP;
  real<lower=0, upper=1> p_DP;
  
  // residual standard deviation
  vector<lower=0>[V] sigma_cause_PC;
  vector<lower=0>[V] sigma_cause_DP;
}

transformed parameters {
  vector<lower=0, upper=1>[V] K_PC;   // causal strength
  vector<lower=0, upper=1>[V] K_DP;

  for (v in 1:V) {
    real p_E = p_PC * (1- (p_PP[v] * (1-p_DP)));
    K_PC[v] = sqrt((p_PC*(1-p_PC))/(p_E*(1-p_E))) * (1-p_PP[v] + p_PP[v]*p_DP);
    K_DP[v] = sqrt((p_DP*(1-p_DP))/(p_E*(1-p_E))) * p_PC * p_PP[v];
  }
}

model {
  // uniform prior over prior probabilities
  p_PC ~ uniform(0, 1);
  p_PP ~ uniform(0, 1);
  p_DP ~ uniform(0, 1);

  // prior over standard deviations
  sigma_cause_PC ~ normal(0, 1);
  sigma_cause_DP ~ normal(0, 1);
  
  if (!prior_only) {
    for (n in 1:N) {
      if (factor[n] == 1)
        cause[n] ~ normal(K_PC[vignette[n]], sigma_cause_PC[vignette[n]]);
      else
        cause[n] ~ normal(K_DP[vignette[n]], sigma_cause_DP[vignette[n]]);
    }
  }
}

generated quantities {
  vector[N] cause_hat;
  vector[N] log_lik;

  for (n in 1:N) {
    if (factor[n] == 1) {
      cause_hat[n] = normal_rng(K_PC[vignette[n]], sigma_cause_PC[vignette[n]]);
      log_lik[n] = normal_lpdf(cause[n] | K_PC[vignette[n]], sigma_cause_PC[vignette[n]]);
    } else {
      cause_hat[n] = normal_rng(K_DP[vignette[n]], sigma_cause_DP[vignette[n]]);
      log_lik[n] = normal_lpdf(cause[n] | K_DP[vignette[n]], sigma_cause_DP[vignette[n]]);
    }
  }
}
