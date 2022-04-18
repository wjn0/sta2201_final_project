library(tidyverse)
library(clock)

extract_state_covid_case_counts <- function(state_sel, start_date, end_date, root_path, remove_unassigned = TRUE) {
  # Extract the COVID case counts for a given selection of states.
  #
  # Args:
  #   state_sel: The list of states (full name, capitalized) to extract.
  #   start_date: The starting date of extraction in day/month/Year format.
  #   end_date: The ending date of extraction in day/month/Year format.
  #   root_path: The path to the `csse_covid_19_daily_reports` data directory.
  #   remove_unassigned: Whether or not to remove cases not assigned to a county.
  #
  # Returns:
  #   A data frame with columns FIPS, county, state, cases, and date.
  country_sel <- c("US")
  start_date <- date_parse(start_date, format = "%d/%m/%Y")
  end_date <- date_parse(end_date, format="%d/%m/%Y")
  
  all_cases <- NULL
  date <- start_date
  while (date <= end_date) {
    path <- paste(root_path, "/", date_format(date, format = "%m-%d-%Y"), ".csv", sep = "")
    
    daily_cases <- read.table(path, sep = ",", quote = '"', header = 1)
    
    daily_cases_filt <- daily_cases %>%
      filter(Country_Region %in% country_sel) %>%
      filter(Province_State %in% state_sel) %>%
      filter(!is.na(FIPS)) %>%
      mutate(county = Admin2, state = Province_State, cases = Confirmed) %>%
      select(FIPS, county, state, cases) %>%
      mutate(date = date)
    if (remove_unassigned) daily_cases_filt <- daily_cases_filt %>% filter(county != "Unassigned")
    
    all_cases <- all_cases %>% bind_rows(daily_cases_filt)
    
    date <- date + 1
  }
  
  return(all_cases)
}

smooth_cumulative_case_counts <- function(cumulative_case_counts) {
  # Given a series of cumulative case counts with potential negative values,
  # smooth it to a series of daily case counts.
  #
  # When we encounter a negative value, we assume the "false reports" that led
  # to the correction are evenly distributed across the prior days, weighted by
  # the number of cases on that day.
  N <- length(cumulative_case_counts)
  if (N < 2) return(NA)
  
  x <- cumulative_case_counts[2:N] - cumulative_case_counts[1:(N - 1)]
  for (i in 1:N) {
    if (!is.na(x[i]) && (x[i] < 0)) {
      to_distribute <- x[i]  # number of negative cases to even out
      weights <- x[1:(i - 1)] / sum(x[1:(i - 1)])
      distribution <- floor(weights * to_distribute)
      x[1:(i - 1)] <- x[1:(i - 1)] + distribution
      x[i] <- 0
    }
  }
  
  return(c(NA, x))
}

extract_smoothed_daily_cases <- function(state_sel, start_date, end_date, root_path, remove_unassigned = TRUE) {
  all_cases <- extract_state_covid_case_counts(state_sel, start_date, end_date, root_path, remove_unassigned)
  smoothed_daily <- all_cases %>%
    group_by(county) %>%
    mutate(cumulative_cases = cases, cases = smooth_cumulative_case_counts(cases))
  
  return(smoothed_daily)
}

extract_smoothed_daily_cases_wide <- function(state_sel, start_date, end_date, root_path, remove_unassigned = TRUE) {
  smoothed_daily <- extract_smoothed_daily_cases(state_sel, start_date, end_date, root_path, remove_unassigned)
  
  smoothed_daily_wide <- smoothed_daily %>%
    pivot_wider(id_cols = date, names_from = county, values_from = cases) %>%
    tail(-1) %>%
    select(where(~!any(is.na(.))))
  
  return(smoothed_daily_wide)
}
