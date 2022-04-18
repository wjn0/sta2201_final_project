data {
  int T; // number of timesteps
  
  real time[T];
  real intervention_time;
  real case_count[T];
}

transformed data {
  real treatment[T];
  
  for(t in 1:T) {
    treatment[t] = (time[t] >= intervention_time);
  }
}

parameters {
  real intercept;
  real time_effect;
  real posttreatment_effect_diff; // non-centered param b/c correlation in posterior hinder test
  
  real<lower=0> epsilon;
}

transformed parameters {
  real mu[T];
  
  for (t in 1:T) {
    mu[t] = (intercept + // global intercept
             (1 - treatment[t]) * time_effect * time[t] + // pre-intervention line
             treatment[t] * (time_effect * intervention_time + (time_effect + posttreatment_effect_diff) * time[t])); // post-intervention line
  }
}

model {
  intercept ~ normal(0, 100);
  time_effect ~ normal(0, 100);
  posttreatment_effect_diff ~ normal(0, 100);
  epsilon ~ normal(0, 100);
  
  case_count ~ normal(mu, epsilon);
}