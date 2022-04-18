# Cacheable model fitting
load_or_do_and_save_fit <- function(fitexpr, filename) {
  full_filename <- paste(filename, ".Rdata", sep = "")
  if (file.exists(full_filename)) {
    fit <- get(load(full_filename))
  }
  else {
    fit <- eval(fitexpr)
    save(fit, file = full_filename)
  }
  
  return(fit)
}
