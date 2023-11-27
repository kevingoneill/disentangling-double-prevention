data {
  int<lower=1> N;                           // number of data points
  int<lower=1> C;                           // number of conditions (2: reversal/adversarial)
  int<lower=1> V;                           // number of vignettes (5)
  int<lower=1> F;                           // number of factors being judged (2: PC/DP)
  array[N] int<lower=1, upper=C> condition; // condition per observation (1: reversal, 2: adversarial)
  array[N] int<lower=1, upper=V> vignette;  // vignette per observation (1: reversal, 2: adversarial)
  array[N] int<lower=1, upper=F> factor;    // factor per observation (1: PC, 2: DP)
  vector<lower=0, upper=1>[N] cause;        // causal judgments
  int<lower=0, upper=1> prior_only;         // sample from just the prior (1) or the posterior (0)?
}

parameters {
  array[V] real<lower=0, upper=1> p_PC;
  array[V] vector<lower=0, upper=1>[C] p_PP;
  array[V] real<lower=0, upper=1> p_DP;

  // shift & scale parameters to align model & participant scales
  real<lower=-1, upper=1> shift;
  real log_scale;
  
  // residual standard deviation
  array[V] vector<lower=0>[C] sigma_cause_PC;
  array[V] vector<lower=0>[C] sigma_cause_DP;
}

transformed parameters {
  real<lower=0> scale = exp(log_scale);
  array[V] vector<lower=0, upper=1>[C] K_PC;   // causal strength
  array[V] vector<lower=0, upper=1>[C] K_DP;

  for (v in 1:V) {
    for (c in 1:C) {
      real p_E = p_PC[v] * (1 - (p_PP[v,c] * (1-p_DP[v])));
      K_PC[v,c] = sqrt((p_PC[v]*(1-p_PC[v]))/(p_E*(1-p_E))) * (1-p_PP[v,c] + p_PP[v,c]*p_DP[v]);
      K_DP[v,c] = sqrt((p_DP[v]*(1-p_DP[v]))/(p_E*(1-p_E))) * p_PC[v] * p_PP[v,c];
    }
  }
}

model {
  for (v in 1:V) {
    // uniform prior over prior probabilities
    p_PC[v] ~ uniform(0, 1);
    p_PP[v] ~ uniform(0, 1);
    p_DP[v] ~ uniform(0, 1);

    // prior over standard deviations
    sigma_cause_PC[v] ~ normal(0, 1);
    sigma_cause_DP[v] ~ normal(0, 1);
  }

  // prior over shift/scale
  shift ~ std_normal();
  log_scale ~ std_normal();  
  
  if (!prior_only) {
    for (n in 1:N) {
      if (factor[n] == 1)
        cause[n] ~ normal(K_PC[vignette[n], condition[n]]*scale + shift,
			  sigma_cause_PC[vignette[n], condition[n]]);
      else
        cause[n] ~ normal(K_DP[vignette[n], condition[n]]*scale + shift,
			  sigma_cause_DP[vignette[n], condition[n]]);
    }
  }
}

generated quantities {
  vector[N] cause_hat;
  vector[N] log_lik;

  for (n in 1:N) {
    if (factor[n] == 1) {
      cause_hat[n] = normal_rng(K_PC[vignette[n], condition[n]]*scale + shift,
				sigma_cause_PC[vignette[n], condition[n]]);
      log_lik[n] = normal_lpdf(cause[n] | K_PC[vignette[n], condition[n]]*scale + shift,
			       sigma_cause_PC[vignette[n], condition[n]]);
    } else {
      cause_hat[n] = normal_rng(K_DP[vignette[n], condition[n]]*scale + shift,
				sigma_cause_DP[vignette[n], condition[n]]);
      log_lik[n] = normal_lpdf(cause[n] | K_DP[vignette[n], condition[n]]*scale + shift,
			       sigma_cause_DP[vignette[n], condition[n]]);
    }
  }
}
