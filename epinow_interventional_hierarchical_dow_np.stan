data {
  int<lower=1> T; // number of timesteps
  int<lower=1> S; // number of series
  int<lower=1> tau; // window size for windowed computations
  int<lower=3> K; // number of spline knots
  int<lower=1> period; // probably 7 days for 1 week
  
  row_vector<lower=0, upper=1>[tau] generation_time; // generation time distribution
  row_vector<lower=0, upper=1>[tau] reporting_delay; // reporting delay distribution
  
  int case_counts[S, T];
  matrix[T - tau, K] B;
  matrix[S, T] treatment;
  int treatment_idx[S];
}

transformed data {
  int omega_idx[T - tau];
  
  for (t in (tau + 1):T) {
    omega_idx[t - tau] = (t % period) + 1;
  }
}

parameters {
  // Time-varying reproduction number R(t).
  real<lower=0> sigma;
  real<lower=0> epsilon;
  vector[K] alpha_mu;
  vector[K] alpha_treat_mu;
  vector[S] alpha_diff;
  vector[S] alpha_treat_diff;
  //matrix[K, S] alpha;
  //matrix[K, S] alpha_treat;
  
  // Treatment effect on reproduction number.
  //real treatment_effect_mean;
  //real<lower=0> treatment_effect_sd;
  //vector[S] treatment_effect;
  
  // Unobserved initial incidence.
  real<lower=0> initial_incidence[S];
  
  // Multiplicative day of the week effect
  real<lower=0> omega[period];
  
  // Overdispersion parameters for negative binomial likelihood
  real<lower=0> phi[S];
}

transformed parameters {
  real incidence[S, T];
  real case_count_mean[S, T - tau];
  matrix[S, T - tau] log_reproduction;
  matrix[S, T - tau] reproduction;
  matrix[K, S] alpha;
  matrix[K, S] alpha_treat;
  
  for(s in 1:S) {
    for (k in 1:K) {
      alpha[k, s] = alpha_mu[k] + alpha_diff[s];
      alpha_treat[k, s] = alpha_treat_mu[k] + alpha_treat_diff[s];
    }
  }
  log_reproduction = ((B * alpha) .* (1 - treatment[, (tau + 1):T]') +
                      (B * alpha_treat) .* treatment[, (tau + 1):T]')';
  reproduction = exp(log_reproduction);
  
  for (t in 1:tau) {
    incidence[, t] = initial_incidence;
  }
  for (s in 1:S) {
    for (t in (tau + 1):T) {
      incidence[s, t] = reproduction[s, t - tau] * (generation_time * to_vector(incidence[s, (t - tau):(t - 1)]));
    }
  }
  
  for (s in 1:S) {
    for (t in (tau + 1):T) {
      case_count_mean[s, t - tau] = reporting_delay * to_vector(incidence[s, (t - tau):(t - 1)]);
    }
  }
}

model {
  // Prior over reproduction number - can be replaced with any functional form
  // we like.
  sigma ~ normal(0, 1);
  alpha_mu[1:2] ~ normal(0, 1);
  alpha_mu[3:K] ~ normal(2 * alpha_mu[2:(K - 1)] - alpha_mu[1:(K - 2)], sigma);
  alpha_treat_mu[1:2] ~ normal(0, 1);
  alpha_treat_mu[3:K] ~ normal(2 * alpha_treat_mu[2:(K - 1)] - alpha_treat_mu[1:(K - 2)], sigma);
  
  epsilon ~ normal(0, 1);
  alpha_diff ~ normal(0, epsilon);
  alpha_treat_diff ~ normal(0, epsilon);
  //for (s in 1:S) {
  //  alpha[, s] ~ normal(alpha_mu, epsilon);
  //  alpha_treat[, s] ~ normal(alpha_treat_mu, epsilon);
  //}
  
  // Treatment effect
  //treatment_effect_mean ~ normal(0, 1);
  //treatment_effect_sd ~ normal(0, 1);
  //treatment_effect ~ normal(treatment_effect_mean, treatment_effect_sd);
  
  // Prior over initial incidence
  initial_incidence ~ lognormal(3, 3);
  
  // Day-of-the week effect
  omega ~ normal(1, 1);
  
  // Overdispersion prior
  phi ~ exponential(1);
  
  // Likelihood
  for (s in 1:S) {
    case_counts[s, (tau + 1):T] ~ neg_binomial_2(to_vector(case_count_mean[s, ]) .* to_vector(omega[omega_idx]) + 1, phi[s]);
  }
}

generated quantities {
  real treatment_effect_mean[S];
  real treatment_effect[S];
  
  for (s in 1:S) {
    treatment_effect_mean[s] = B[treatment_idx[s] - tau, ] * (alpha_treat_mu - alpha_mu);
    treatment_effect[s] = B[treatment_idx[s] - tau, ] * (alpha_treat[, s] - alpha[, s]);
  }
}
