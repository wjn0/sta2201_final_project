library(rstan)
library(zoo)

source("simulation.R")
source("getsplines.R")
source("modeling.R")

# Simulate the data that we use for the initial model tests.
set.seed(10)
intervention_time <- 100
simulated <- simulate(T = 365, discont = 100)
save(simulated, file = "simulated.Rdata")
tau <- length(simulated$w)

# RD Gaussian fit
rd_gaussian_data <- list(T = length(simulated$C),
                         time = seq(1, length(simulated$C)),
                         intervention_time = 100,
                         case_count = simulated$C)
rd_gaussian_fit <- load_or_do_and_save_fit(
  quote(stan(file = "rd_gaussian.stan", data = rd_gaussian_data)),
  "rd_gaussian"
)

# RD Gaussian small window
window_size <- 30
time_idx <- seq(intervention_time - window_size, intervention_time + window_size - 1)
rd_gaussian_data_small_window <- list(T = window_size * 2,
                                      time = time_idx,
                                      intervention_time = intervention_time,
                                      case_count = simulated$C[time_idx])
rd_gaussian_small_window_fit <- load_or_do_and_save_fit(
  quote(stan(file = "rd_gaussian.stan", data = rd_gaussian_data_small_window)),
  "rd_gaussian_small_window"
)

# RK Gaussian fit - still using small window
window_size <- 30
time_idx <- seq(intervention_time - window_size, intervention_time + window_size - 1)
rk_gaussian_data_small_window <- list(T = window_size * 2,
                                      time = time_idx,
                                      intervention_time = intervention_time,
                                      case_count = rollmean(simulated$C, 7, na.pad = TRUE)[time_idx])
rk_gaussian_fit <- load_or_do_and_save_fit(
  quote(stan(file = "rk_gaussian.stan", data = rd_gaussian_data_small_window)),
  "rk_gaussian"
)

# RK Gaussian w/ effect in wrong place
rk_gaussian_misplaced_data <- list(T = window_size * 2,
                                   time = time_idx,
                                   intervention_time = intervention_time + 14,
                                   case_count = rollmean(simulated$C, 7, na.pad = TRUE)[time_idx])
rk_gaussian_misplaced_fit <- load_or_do_and_save_fit(
  quote(stan(file = "rk_gaussian.stan", data = rk_gaussian_misplaced_data)),
  "rk_gaussian_misplaced"
)

# EpiNow - interventional
tau <- length(simulated$w)
T <- length(simulated$C)
B <- getsplines(1:length(simulated$C), I = 28)$B.ik[(tau + 1):T, ]
epinow_data <- list(T = T,
                    tau = tau,
                    omega = c(rep(c(simulated$omega[2:7], simulated$omega[1]), 52), 0),
                    generation_time = rev(simulated$w),  # stan does not allow reverse indexing
                    reporting_delay = rev(simulated$xi),  # stan does not allow reverse indexing
                    case_counts = simulated$C,
                    K = ncol(B),
                    B = B,
                    treatment = (seq(1:length(simulated$C)) >= 100))

epinow_interventional_fit <- load_or_do_and_save_fit(
  quote(stan(file = "epinow_interventional.stan", data = epinow_data, cores = 4)),
  "epinow_interventional"
)

# EpiNow - hierarchical
set.seed(6)
discont_loc <- 60
simulated_h <- simulate_hierarchical(Q = 10, T = 120, discont = discont_loc)
save(simulated_h, file = "simulated_hierarchical.Rdata")

case_counts <- matrix(unlist(simulated_h$Cs), nrow = length(simulated_h$Cs), byrow = TRUE)

tau <- length(simulated_h$w)
T <- ncol(case_counts)
S <- nrow(case_counts)

B <- getsplines((tau + 1):T, I = 10)$B.ik

omega <- c(rep(c(simulated_h$omega[2:7], simulated_h$omega[1]), 52), 0)

treatment <- matrix(rep((seq(1:length(simulated_h$Cs[[1]])) >= discont_loc), S),
                    nrow = length(simulated_h$Cs),
                    byrow = TRUE) * 1

epinow_h_data <- list(S = S,
                      T = T,
                      tau = tau,
                      generation_time = rev(simulated_h$w),  # stan does not allow reverse indexing
                      reporting_delay = rev(simulated_h$xi),  # stan does not allow reverse indexing
                      case_counts = case_counts,
                      K = ncol(B),
                      B = B,
                      treatment = treatment,
                      period = 7,
                      treatment_idx = rep(discont_loc, S))
save(epinow_h_data, file = "epinow_h_data.Rdata")

epinow_interventional_h_fit <- load_or_do_and_save_fit(
  quote(stan(file = "epinow_interventional_hierarchical_dow_np.stan", data = epinow_h_data, cores = 4)),
  "epinow_interventional_h"
)
