/*
 * Modeling of cover class data
 * using regularized incomplete beta function
 */

functions {
 /*
  * Return the log probability that the cover class is oberved under the
  * given paramters a and b.
  *
  * @param Y  Observed class
  * @param CP Cut points
  * @param a  Parameter of beta distribution
  * @param b  Parameter of beta distribution
  */
  real coverclass_lpmf(int Y, vector CP, real a, real b) {
    int n_cls = num_elements(CP) + 1;
    real gamma;

    if (Y <= 1) {  // 0 or 1
      gamma =  inc_beta(a, b, CP[1]);
    } else if(Y >= 2 && Y < n_cls) {
      gamma = inc_beta(a, b, CP[Y])
              - inc_beta(a, b, CP[Y - 1]);
    } else {
      gamma = 1 - inc_beta(a, b, CP[n_cls - 1]);
    }
    return bernoulli_lpmf(1 | gamma);
  }
  
 /*
  * Return cover class randomly given cutpoints and paramters a and b.
  *
  * @param CP     Cut points
  * @param n_cls  Number of classes
  * @param a      Parameter of beta distribution
  * @param b      Parameter of beta distribution
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
  int<lower = 1> N_cls;                       // Number of classes
  int<lower = 1> N_t;                         // Number of observations
  int<lower = 0, upper = N_cls> Y[N_t];       // Observed cover class
  vector<lower = 0, upper = 1>[N_cls - 1] CP; // Cut points
}

parameters {
  real<lower = 0, upper = 1> delta;       // Intra-quad corr.
                                          //   or uncertainty
  vector[N_t] theta_raw;                  // Reparameterized state changes
  real<lower = 0> sigma;                  // S.D. of temporal changes
}

transformed parameters {
  vector[N_t] theta;                      // Latent state
  vector[N_t] mu;                         // Mean proportion of cover

  // System model (reparameterized)
  theta[1:2] = theta_raw[1:2];
  for (t in 3:N_t)
    theta[t] = 2 * theta[t - 1] - theta[t - 2] + theta_raw[t] * sigma;
  mu = inv_logit(theta);
}

model {
  // System model (reparameterized)
  theta_raw[1:2] ~ normal(0, 2.5);    // Prior
  theta_raw[3:N_t] ~ std_normal();
  
  // Observation
  for (t in 1:N_t) {
    real a = mu[t] / delta - mu[t];
    real b = (1 - mu[t]) * (1 - delta) / delta;

    Y[t] ~ coverclass(CP, a, b);
  }

  // Prior
  sigma ~ normal(0, 2.5);
}

generated quantities {
  int yrep[N_t]; // for posterior predictive check
  
  for (t in 1:N_t) {
    real a = mu[t] / delta - mu[t];
    real b = (1 - mu[t]) * (1 - delta) / delta;

    yrep[t] = coverclass_rng(CP, N_cls, a, b);
  }
}
