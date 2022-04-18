library(rstan)
library(here)
library(clock)

source("simulation.R")
source("data.R")
source("getsplines.R")
source("modeling.R")

# ensure no functions are masked /facepalm R why are you like this
library(tidyverse)

# options for speed
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Load data
root_path <- here("data/jhu_covid/COVID-19-master/csse_covid_19_data/csse_covid_19_daily_reports")
cases <- extract_smoothed_daily_cases_wide("California", "01/01/2022", "15/04/2022", root_path)
cases_sub <- cases[
    , c("date", "Kern", "Fresno", "Alameda", "Sacramento", "Santa Clara", "San Bernardino",
        "Riverside", "Orange", "San Diego", "Los Angeles")
]

# Structure model inputs
w <- make_discrete_lognormal(mean = 3.6, sd = 3.1)  # from EpiNow paper - kinda
xi <- make_discrete_lognormal(mean = 5.2, sd = 1.52)  # from EpiNow paper
case_counts <- t(as.matrix(cases_sub[, 2:ncol(cases_sub)]))
tau <- length(w)
T <- ncol(case_counts)
S <- nrow(case_counts)
B <- getsplines((tau + 1):T, I = 10)$B.ik

# The date of lifting the mandate
treatment_vec <- cases_sub$date >= date_parse("01-03-2022", format = "%d-%m-%Y")
treatment <- matrix(rep(treatment_vec, S), nrow = S, byrow = TRUE) * 1

epinow_data <- list(S = S,
                    T = T,
                    tau = tau,
                    generation_time = rev(w),  # stan does not allow reverse indexing
                    reporting_delay = rev(xi),  # stan does not allow reverse indexing
                    case_counts = case_counts,
                    K = ncol(B),
                    B = B,
                    treatment = treatment,
                    period = 7,
                    treatment_idx = rep(which.max(treatment_vec), S))

california_fit <- load_or_do_and_save_fit(
    quote(stan(file = "epinow_interventional_hierarchical_dow_np.stan",
               data = epinow_data)),
    "california_fit"
)
