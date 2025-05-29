# Function to fit natural spline-based linear models to null and alternate
# hypotheses H0 and H1
fit_hypothesis_splines <- function(n, x_tbl, d_f_total, silent = FALSE){
  d_f <- d_f_total[n]
  num_total <- length(d_f_total)
  # Build linear model formula
  spline_formula <- paste0("quantity ~ splines::ns(Temp, df = ", d_f, ")")
  alt_spline_formula <- paste0(spline_formula, " * Condition")

  if(!silent) {
    progress <-  floor(5 * (n-1) / num_total) * 2
    cat("\r|", strrep("=", progress), strrep(" ", 10 - progress), "| ",
        d_f, " degrees of freedom: H-null...", sep = "")
  }
  x_tbl <- split(x_tbl, ~ Protein_ID)
  null_fit <- get_fit_details(x_tbl, spline_formula, "H-null")

  if(!silent) {
    progress <-  floor(5 * (n-1) / num_total) * 2 + 1
    cat("\r|", strrep("=", progress), strrep(" ", 10 - progress), "| ",
        d_f, " degrees of freedom: H-alt... ", sep = "")
  }
  alt_fit <- get_fit_details(x_tbl, alt_spline_formula, "H-alt")

  if(!silent & n == num_total) {
    cat("\r|", strrep("=", 10), "| ",
        "Done!                               \n", sep = "")
  }

  combined_fit <- rbind(null_fit, alt_fit)
  combined_fit$degrees_of_freedom <- d_f
  combined_fit
}


# Transform model into dataframe of details
get_spline_model_details <- function(sp_fit){
  if (inherits(sp_fit, "lm")){

    # Residual sum of squares
    fit_RSS <- stats::deviance(sp_fit)

    # Residual standard deviation (sigma)
    fit_sigma <- stats::sigma(sp_fit)

    # Number of observations / coefficients
    fit_n_obs <- stats::nobs(sp_fit)
    fit_n_coeffs <- length(stats::coef(sp_fit))

    # Log-likelihood / Corrected Akaike information criterion (AICc)
    fit_log_lik <- stats::logLik(sp_fit)

    k <- attr(fit_log_lik, "df")
    fit_AICc <- stats::AIC(sp_fit) + 2 * k * (k + 1)/(fit_n_obs - k - 1)

    fit_log_lik <- as.numeric(fit_log_lik)


  } else {

    fit_RSS <- NA
    fit_sigma <- NA
    fit_n_coeffs <- NA
    fit_n_obs <- NA
    fit_log_lik <- NA
    fit_AICc <- NA

  }

  fit_details <-
    data.frame(RSS = fit_RSS,
               sigma = fit_sigma,
               n_coeffs = fit_n_coeffs,
               n_obs = fit_n_obs,
               log_lik = fit_log_lik,
               AICc = fit_AICc)
}

# Function to get and decorate model details per protein
get_fit_details <- function(x_tbl, spline_formula, hypothesis){
  fit_details <-
    lapply(x_tbl, get_natural_spline_lm, formula = spline_formula) |>
    Reduce(rbind, x = _)
  fit_details$success <- !is.na(fit_details$RSS)
  fit_details$hypothesis <- hypothesis

  fit_details
}

# Function to fit linear model natural spline
get_natural_spline_lm <- function(x_tbl, formula, silent = TRUE){
  spline_lm <- try(stats::lm(formula, x_tbl), silent = silent)
  if(inherits(spline_lm, "try-error")) spline_lm <- NULL
  lm_details <- get_spline_model_details(spline_lm)

  cbind(data.frame(Protein_ID = x_tbl$Protein_ID[1]), lm_details)
}
