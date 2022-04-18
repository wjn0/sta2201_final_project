library(MASS)

# Functions for simulating data (roughly) according to the EpiNow 2 model described
# in https://wellcomeopenresearch.org/articles/5-112
#
# We deviate in the simulation of the R(t), using a different structure for the GP
# (I think).

make_discrete_lognormal <- function(mean, sd, tau = 10) {
  # Make a discrete distribution from a log-normal. Used for incubation period and
  # generation time distributions.
  mu <- log(mean^2 / sqrt(mean^2 + sd^2))
  sigma <- sqrt(log(1 + sd^2 / mean^2))
  
  cdf <- c()
  discretization <- c()
  for (t in 1:tau) {
    if (t == tau) cdf[t] <- 1
    else cdf[t] <- pnorm(log(t), mean = mu, sd = sigma)
  
    if (t > 1) discretization[t] <- cdf[t] - cdf[t - 1]
    else discretization[t] <- cdf[t]
  }
  
  return(rev(discretization))
}

simulate_R <- function(T = 100, mean = 0.3, sd = 0.4, discont = FALSE, discont_magnitude = -0.2) {
  # Simulates R(t) from a zero-mean Gaussian process with a squared-exponential
  # kernel, and transforms via `exp` to obtain observations on $[0, \infty]$.
  #
  # Args:
  #     T: The number of days to simulate.
  #     sd: The standard deviation of the Gaussian process.
  #     discont: If FALSE, no discontinuity is simulated. If a number, must be
  #              less than T and greater than 1 and be the day of the discontinuity.
  #
  # Returns:
  #     R: R(t), a vector of length T.
  
  # Simulate from GP
  ts <- seq(1, T)
  D <- abs(outer(ts, ts, FUN = "-"))  # euclidean distance between times
  sigma <- median(D)  # kernel bandwidth by median heuristic
  K <- exp(-D**2 / (2*sigma**2))
  log_R <- mvrnorm(mu = rep(0, T), Sigma = K) * sd + mean
  
  # Simulate discontinuity
  if (discont) {
    stopifnot(discont < T, discont > 1)
    
    log_R[discont:T] <- log_R[discont:T] + discont_magnitude
  }
  
  # Transform to $[0, \infty]$
  R <- exp(log_R)
  
  return(R)
}

simulate_w <- function(tau = 14) {
  # Simulate the probability mass function for transmission. This is the
  # probability at an (unobserved) infection on day t is a transmission event
  # from an (unobserved) infection spawned on day t - tau.
  #
  # We opt for a simple unimodal mass function centred at tau / 2.
  #
  # Args:
  #     tau: The window for transmission.
  #
  # Returns:
  #     w: The probability mass function for transmission, a vector of length
  #        tau.
  xs <- seq(-3, 3, by = 6 / (tau - 1))
  w_unnorm <- dnorm(xs)
  w <- w_unnorm / sum(w_unnorm)
  
  return(w)
}

simulate_I <- function(r, w, initial_mean = 10) {
  # Simulate the latent infections based on the reproduction number and
  # transmission probabilities.
  #
  # Args:
  #     r: The time-varying reproduction number, of length T.
  #     w: The transmission probability mass function, of length tau.
  #     initial_mean: The infections on the first tau days will be drawn from a
  #                   Poisson with this mean.
  #
  # Returns:
  #     I: The latent infections, of length T + tau.
  T <- length(r)
  tau <- length(w)
  
  # Simulate initial infections ~ Poisson(initial_mean)
  I <- rep(0, T + tau)
  I[1:tau] <- rpois(tau, initial_mean)
  
  # Simulate the rest of the infections
  for (t in 1:T) {
    I[tau + t] <- r[t] * (w %*% I[(t + tau - 1):t])
  }
  
  return(I)
}

simulate_xi <- function(tau = 14) {
  # Simulate reporting delay/incubation period mass function. We adopt the same
  # technique as we did for w.
  #
  # Args:
  #     tau: The window size.
  #
  # Returns:
  #     xi: The reporting delay/incubation period mass function, of length tau.
  xi <- simulate_w(tau = tau)
  
  return(xi)
}

simulate_D <- function(xi, I) {
  tau <- length(xi)
  T <- length(I) - tau
  
  D <- rep(0, T)
  for (t in 1:T) {
    D[t] <- xi %*% I[(t + tau - 1):t]
  }
  
  return(D)
}

simulate_C <- function(D, omega, phi) {
  stopifnot(length(omega) == 7)
  T <- length(D)
  
  omega_t <- omega[(seq(1, T) %% 7) + 1]
  mean <- D * omega_t
  C <- rnegbin(n = T, mean, phi)
  
  return(C)
}

simulate_omega <- function() {
  omega <- c(rep(0, 2), 3, runif(n = 4) / 2 + 0.5)
  
  return(omega)
}

simulate <- function(T = 100, discont = FALSE) {
  R <- simulate_R(T = T, discont = discont)
  w <- simulate_w()
  I <- simulate_I(R, w)
  xi <- simulate_xi()
  D <- simulate_D(xi, I)
  omega <- simulate_omega()
  phi <- rexp(n = 1, rate = 1)
  C <- simulate_C(D, omega, phi)
  
  ret <- list(R = R, w = w, I = I, xi = xi, D = D, omega = omega, phi = phi, C = C)
  return(ret)
}

simulate_hierarchical <- function(Q = 15, T = 365, discont = 100, discont_magnitude_mean = -0.2, discont_magnitude_sd = 0.05) {
  # Simulate base R(t)
  base_R <- simulate_R(T = T, mean = 0.3, discont = FALSE)
  discont_magnitudes <- rnorm(n = Q, mean = discont_magnitude_mean, sd = discont_magnitude_sd)
  
  # Simulate shared parameters
  w <- make_discrete_lognormal(mean = 3.6, sd = 3.1)  # from EpiNow paper - kinda
  xi <- make_discrete_lognormal(mean = 5.2, sd = 1.52)  # from EpiNow paper
  omega <- simulate_omega()
  
  # Simulate specific parameters
  phi <- rexp(n = Q, rate = 1)
  Rs <- list()
  Is <- list()
  Ds <- list()
  Cs <- list()
  for (q in 1:Q) {
    Rs[[q]] <- exp(log(base_R) + discont_magnitudes[q] * (seq(1:T) >= discont))
    Is[[q]] <- simulate_I(Rs[[q]], w)
    Ds[[q]] <- simulate_D(xi, Is[[q]])
    Cs[[q]] <- simulate_C(Ds[[q]], omega, phi[q])
  }
  
  ret <- list(w = w, xi = xi, base_R = base_R, phi = phi, Rs = Rs, Is = Is, omega = omega,
              Ds = Ds, Cs = Cs, discont_magnitudes = discont_magnitudes)
  return(ret)
}
