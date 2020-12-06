/*
 * A model to estimate plant cover
 * from cover class data
 * using regularized incomplete beta function
 */

functions {
 /*
  * Return the log probability the cover class Y is observed
  * given the cut points CP and thebeta distribution parameters
  * a and b
  *
  * @param Y  Observed cover class (int)
  * @param CP Cut points (vector)
  * @param a  Parameteer of the beta distribution (real)
  * @param b  Parameteer of the beta distribution (real)
  *
  * @return Log probability that Y is observed
  */
  real coverclass_lpmf(int Y, vector CP, real a, real b) {
    int n_class;
    real gamma;

    n_class = num_elements(CP) + 1;
    if (Y <= 1) {  // 0 or 1
      gamma =  inc_beta(a, b, CP[1]);
    } else if(Y >= 2 && Y < n_class) {
      gamma = inc_beta(a, b, CP[Y])
              - inc_beta(a, b, CP[Y - 1]);
    } else {
      gamma = 1 - inc_beta(a, b, CP[n_class - 1]);
    }
    return bernoulli_lpmf(1 | gamma);
  }

 /*
  * Return the cover class randomly given the cut points CP,
  * number of classes, and the beta distribution parameters
  * a and b
  *
  * @param CP    Cut points (vector)
  * @param n_cls Number of classes (int)
  * @param a     Parameteer of the beta distribution (real)
  * @param b     Parameteer of the beta distribution (real)
  *
  * @return Log probability that Y is observed
  */
  int coverclass_rng(vector CP, int n_cls, real a, real b) {
    vector[n_cls] pr;
    int y;
    
    pr[1] = inc_beta(a, b, CP[1]);
    for (i in 2:(n_cls - 1))
      pr[i] = inc_beta(a, b, CP[i]) - inc_beta(a, b, CP[i - 1]);
    pr[n_cls] = 1 - inc_beta(a, b, CP[n_cls - 1]);
    y = categorical_rng(pr);
    return y;
  }
}

data {
  int<lower = 1> N_q;                           // Number of quadrats
  int<lower = 1> N_y;                           // Number of years
  int<lower = 1> N_cls;                         // Number of classes
  int<lower = 1> N_obs;                         // Number of observed years
  int<lower = 1> Obs_y[N_obs];                  // Observed years
  vector<lower = 0, upper = 1>[N_cls - 1] CP;   // Cut points
  int<lower = 0, upper = N_cls> Y[N_obs, N_q];  // Cover class data
}

parameters {
  matrix[N_obs, N_q] r_raw;          // Spatial random effect (reparameterized)
  vector[N_y] theta_raw;             // Latent state (reparameterized)
  real<lower = 0, upper = 1> delta;  // Uncertainty or intra-quadrat corr.
  real<lower = 0> sigma[2];          // Standard deviations of
                                     //   spatial random effect: sigma[1]
                                     //   and temporal variation: sigma[2]
}

transformed parameters {
  matrix[N_obs, N_q] r;              // Spatial random effect
  vector[N_y] theta;                 // Latent state
  
  // System model (reparameterized)
  // Spatial random effect
  for (i in 1:N_obs) {
    r[i, 1] = r_raw[i, 1];
    for (q in 2:N_q) {
      r[i, q] = r[i, q - 1] + r_raw[i, q] * sigma[1];
    }
    r[i] = r[i] - mean(r[i]);
  }

  // Logit of cover proportion
  theta[1:2] = theta_raw[1:2];
  for (t in 3:N_y)
    theta[t] = 2 * theta[t - 1] - theta[t - 2] + theta_raw[t] * sigma[2];
}

model {
  // Variation in spatial random effect
  for (i in 1:N_obs) {
    r_raw[i, 1] ~ normal(0, 2.5);     // Prior
    r_raw[i, 2:N_q] ~ std_normal();
  }

  // Variation in system model
  theta_raw[1:2] ~ normal(0, 2.5);    // Prior
  theta_raw[3:N_y] ~ std_normal();

  // Observation model
  for (i in 1:N_obs) {
    int y = Obs_y[i];
    for (q in 1:N_q) {
        real mu = inv_logit(theta[y] + r[i, q]);
        real alpha = mu / delta - mu;
        real beta = (1 - mu) * (1 - delta) / delta;

        Y[i, q] ~ coverclass(CP, alpha, beta);
    }
  }

  // Weakly informative priors
  sigma ~ normal(0, 2.5);
}

generated quantities {
  vector[N_y] phi = inv_logit(theta);
  int yrep[N_obs, N_q];

  for (i in 1:N_obs) {
    int y = Obs_y[i];
    for (q in 1:N_q) {
      real mu = inv_logit(theta[y] + r[i, q]);
      real alpha = mu / delta - mu;
      real beta = (1 - mu) * (1 - delta) / delta;

      yrep[i, q] = coverclass_rng(CP, N_cls, alpha, beta);
    }
  }
}
