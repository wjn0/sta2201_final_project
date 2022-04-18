data {
  int<lower=1> T; // number of timesteps
  int<lower=1> tau; // window size for windowed computations
  int<lower=3> K; // number of spline knots
  
  real<lower=0> omega[T]; // multiplicative day of the week effect
  
  row_vector<lower=0, upper=1>[tau] generation_time; // generation time distribution
  row_vector<lower=0, upper=1>[tau] reporting_delay; // reporting delay distribution
  
  int case_counts[T];
  matrix[T - tau, K] B;
  vector[T] treatment;
}

parameters {
  // Time-varying reproduction number R(t).
  real<lower=0> sigma;
  vector[K] alpha;
  
  // Treatment effect on reproduction number.
  real treatment_effect;
  
  // Unobserved initial incidence.
  real<lower=0> initial_incidence;
  
  // Overdispersion parameters for negative binomial likelihood
  real<lower=0> phi;
}

transformed parameters {
  real incidence[T];
  real case_count_mean[T - tau];
  vector[T - tau] log_reproduction = B * alpha + treatment_effect * treatment[(tau + 1):T];
  vector[T - tau] reproduction = exp(log_reproduction);
  
  for (t in 1:tau) {
    incidence[t] = initial_incidence;
  }
  for (t in (tau + 1):T) {
    incidence[t] = reproduction[t - tau] * (generation_time * to_vector(incidence[(t - tau):(t - 1)]));
  }
  
  for (t in (tau + 1):T) {
    case_count_mean[t - tau] = reporting_delay * to_vector(incidence[(t - tau):(t - 1)]);
  }
}

model {
  // Prior over reproduction number - can be replaced with any functional form
  // we like.
  sigma ~ normal(0, 1);
  alpha[1:2] ~ normal(0, 1);
  alpha[3:K] ~ normal(2 * alpha[2:(K - 1)] - alpha[1:(K - 2)], sigma);
  
  // TODO this may break things... did not test but got reasonable results w/o it!?
  treatment_effect ~ normal(0, 1);
  
  // Prior over initial incidence
  initial_incidence ~ lognormal(3, 3);
  
  // Overdispersion prior
  phi ~ exponential(1);
  
  // Likelihood
  case_counts[(tau + 1):T] ~ neg_binomial_2(to_vector(case_count_mean) .* to_vector(omega[(tau + 1):T]) + 1, phi);
}