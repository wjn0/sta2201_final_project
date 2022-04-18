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
  real treatment_effect;
  real time_effect;
  
  real<lower=0> epsilon;
}

transformed parameters {
  real mu[T];
  
  for (t in 1:T) {
    mu[t] = intercept + time_effect * time[t] + treatment_effect * treatment[t];
  }
}

model {
  intercept ~ normal(0, 100);
  treatment_effect ~ normal(0, 100);
  time_effect ~ normal(0, 100);
  epsilon ~ normal(0, 100);
  
  case_count ~ normal(mu, epsilon);
}