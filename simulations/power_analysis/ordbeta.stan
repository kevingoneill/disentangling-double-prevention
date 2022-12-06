functions {
  real ord_beta_reg_lpdf(real y, real mu, real phi, vector cutpoints) {  
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
  int<lower=1> N;       // number of data points
  int<lower=1> K;       // number of predictors
  matrix[N, K] X;       // design matrix
  vector[N] y;          // response variable
}

parameters {
  // parameters for ordered beta
  ordered[2] cutpoints;           // bounds to force 0/1 values

  // coefficients
  vector[K] b;
  vector[K] b_phi;
}

model {
  // priors for intercepts and coefficients
  b[1] ~ student_t(3, 0, 2.5);
  b_phi[1] ~ student_t(3, 0, 2.5);
  for (k in 2:K) {
    b ~ normal(0, 1);
    b_phi ~ normal(0, 1);
  }
  
  cutpoints ~ normal(0, 10);
  
  if (!prior_only) {
    for (n in 1:N)
      target += ord_beta_reg_lpdf(y[n] | inv_logit(X[n,]*b), exp(X[n,]*b_phi), cutpoints);
  }
}
