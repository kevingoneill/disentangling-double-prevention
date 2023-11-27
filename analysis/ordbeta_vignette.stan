functions {
  /*
   * log probability density function of ordered beta distribution
   * Args:
   *   y: a single observation in [0, 1]
   *   mu: the mean of the beta component
   *   phi: the precision of the beta component
   *   cutpoints: two cutpoints on the latent distribution to model 0's and 1's
   */
  real ord_beta_lpdf(real y, real mu, real phi, vector cutpoints) {  
    //auxiliary variables
    real mu_logit = logit(mu);
        
    if (y == 0) {
      return log1m_inv_logit(mu_logit - cutpoints[1]);
    } else if (y == 1) {
      return log_inv_logit(mu_logit - cutpoints[2]);
    } else {
      return log(inv_logit(mu_logit - cutpoints[1]) - inv_logit(mu_logit - cutpoints[2])) +
	beta_proportion_lpdf(y | mu, phi);
    }
  }
}

data {
  int<lower=0, upper=1> prior_only; // sample from the prior?
  int<lower=1> N;                   // number of data points
  int<lower=1> K;                   // number of predictors
  matrix[N, K] X;                   // design matrix
  int<lower=1> V;                   // number of vignettes
  int<lower=1> K_V;                 // number of vignette-level predictors
  array[N] int<lower=1> v;          // vignette indicator
  matrix[N, K_V] X_V;               // vignette-level design matrix
  vector[N] y;                      // response variable

  int<lower=1> N_pred;           // number of predictions
  matrix[N_pred, K] X_pred;      // design matrix for predictions
  matrix[N_pred, K] X_V_pred;    // vignette-level design matrix for predictions
}

transformed data {
  // double the number of vignette-level predictors
  // (one set for mu, one set for phi)
  int Kv = K_V * 2;
}

parameters {
  // parameters for ordered beta
  ordered[2] cutpoints;           // bounds to force 0/1 values

  // population-level coefficients
  vector[K] b;
  vector[K] b_phi;

  // vignette-level effects
  vector<lower=0>[Kv] sd_v;       // vignette-level standard deviations
  matrix[Kv, V] z_v;              // standardized vignette-level effects
  cholesky_factor_corr[Kv] L_v;   // vignette-level effect correlations
}

transformed parameters {
  // actual vignette-level effects
  array[V] vector[K_V] r_v;
  array[V] vector[K_V] r_v_phi;
  
  {
    // pre-index vignette-level effects for efficiency
    matrix[Kv, V] R = diag_pre_multiply(sd_v, L_v) * z_v;
    for (i in 1:V) {
      r_v[i] = R[1:K_V, i];
      r_v_phi[i] = R[(K_V+1):Kv, i];
    }
  }
  
}

model {
  // priors for intercepts and coefficients
  b[1] ~ student_t(3, 0, 2.5);
  b_phi[1] ~ student_t(3, 0, 2.5);
  for (k in 2:K) {
    b ~ normal(0, 1);
    b_phi ~ normal(0, 1);
  }

  // priors for vignette-level effects
  sd_v ~ std_normal();
  to_vector(z_v) ~ std_normal();
  L_v ~ lkj_corr_cholesky(2);
  
  cutpoints ~ normal(0, 1);
  
  if (!prior_only) {
    for (n in 1:N) {
      y[n] ~ ord_beta(inv_logit(X[n,]*b + X_V[n,]*r_v[v[n]]),
		      exp(X[n,]*b_phi + X_V[n,]*r_v_phi[v[n]]),
		      cutpoints);
    }
  }
}

generated quantities {
  // vignette-level effect correlations
  matrix[Kv, Kv] Omega = L_v * transpose(L_v);

  // posterior predictions
  vector[N] y_pred;
  {
    for (n in 1:N) {
      real mu = inv_logit(X[n,]*b + X_V[n,]*r_v[v[n]]);
      real phi = exp(X[n,]*b_phi + X_V[n,]*r_v_phi[v[n]]);
      real p_0 = 1 - inv_logit(logit(mu) - cutpoints[1]);
      real p_01 = inv_logit(logit(mu) - cutpoints[1]) - inv_logit(logit(mu) - cutpoints[2]);
      real p_1 = inv_logit(logit(mu) - cutpoints[2]);
      
      int mixture = categorical_rng([p_0, p_01, p_1]');
      if (mixture == 1)
        y_pred[n] = 0;
      else if (mixture == 2)
        y_pred[n] = beta_proportion_rng(mu, phi);
      else
        y_pred[n] = 1;
    }
  }

  // vignette-average conditional parameters
  vector[N_pred] mu_logit = X_pred * b;
  vector[N_pred] phi_log = X_pred * b_phi;
  vector[N_pred] mu = inv_logit(mu_logit);
  vector[N_pred] phi = exp(phi_log);
  vector[N_pred] p_0 = 1 - inv_logit(logit(mu) - cutpoints[1]);
  vector[N_pred] p_01 = inv_logit(logit(mu) - cutpoints[1]) - inv_logit(logit(mu) - cutpoints[2]);
  vector[N_pred] p_1 = inv_logit(logit(mu) - cutpoints[2]);
  vector[N_pred] e_pred = p_01.*mu + p_1;

  // conditional parameters per vignette
  matrix[N_pred, V] mu_v_logit;
  matrix[N_pred, V] phi_v_log;
  for (i in 1:V) {
    mu_v_logit[, i] = X_pred*b + X_V_pred*r_v[i];
    phi_v_log[, i] = X_pred*b_phi + X_V_pred*r_v_phi[i];
  }

  matrix[N_pred, V] mu_v = inv_logit(mu_v_logit);
  matrix[N_pred, V] phi_v = exp(phi_v_log);
  
  
  matrix[N_pred, V] p_0_v = 1 - inv_logit(logit(mu_v) - cutpoints[1]);
  matrix[N_pred, V] p_01_v = inv_logit(logit(mu_v) - cutpoints[1]) - inv_logit(logit(mu_v) - cutpoints[2]);
  matrix[N_pred, V] p_1_v = inv_logit(logit(mu_v) - cutpoints[2]);
  matrix[N_pred, V] e_pred_v = p_01_v.*mu_v + p_1_v;  
}
