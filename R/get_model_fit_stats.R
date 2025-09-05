# Transform model into dataframe of details
get_model_fit_stats <- function(sp_fit){
  if (inherits(sp_fit, c("lm", "nls"))){

    # Assume sigmoid fit for non-linear model
    if(inherits(sp_fit, "nls")){

      # Get original values
      fit_env <- sp_fit$m$getEnv()
      pre_model_data <- data.frame(x = fit_env$x, y = fit_env$y)

      # R2
      fit_SS_total <-
        sum((pre_model_data$y - mean(pre_model_data$y, na.rm=TRUE))^2, na.rm=TRUE)
      fit_SS_resid <-
        sum((pre_model_data$y - stats::predict(sp_fit))^2 , na.rm=TRUE)
      fit_R_sq <- 1 - fit_SS_resid / fit_SS_total

      # Tm
      fit_pars <- sp_fit$m$getPars()
      a <- fit_pars["a"]
      b <- fit_pars["b"]
      pl <- fit_pars["pl"]
      fit_Tm <-
        suppressWarnings(a / (b - log((1 - pl) / (0.5 - pl) - 1)))
      if(is.na(fit_Tm)) fit_Tm <- NA

      # Tinfl
      fit_Tinterval <- c(min(pre_model_data$x), max(pre_model_data$x))
      fit_Tinfl <-
        get_sigmoid_formula_root(get_sigmoid_TPPTR_formula(2),
                                 a = a, b = b, pl = pl,
                                 x_interval = fit_Tinterval)
      # Slope
      x <- fit_Tinfl
      if(!is.na(fit_Tinfl)){
        fit_slope <- eval(parse(text = get_sigmoid_TPPTR_formula(1)))
      } else {
        fit_slope <- NA
      }

    # Assume spline fit for linear model:
    # no a, b, plateau or sigmoid-derived stats
    } else {
      a <- NA
      b <- NA
      pl <- NA
      fit_Tm <- NA
      fit_Tinfl <- NA
      fit_slope <- NA
      fit_R_sq <- NA
    }

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

    a <- NA
    b <- NA
    pl <- NA
    fit_Tm <- NA
    fit_Tinfl <- NA
    fit_slope <- NA
    fit_R_sq <- NA
    fit_RSS <- NA
    fit_sigma <- NA
    fit_n_coeffs <- NA
    fit_n_obs <- NA
    fit_log_lik <- NA
    fit_AICc <- NA

  }

  fit_details <-
    data.frame(a,
               b,
               plateau = pl,
               melt_point = fit_Tm,
               infl_point = fit_Tinfl,
               slope = fit_slope,
               R_sq = fit_R_sq,
               RSS = fit_RSS,
               sigma = fit_sigma,
               n_coeffs = fit_n_coeffs,
               n_obs = fit_n_obs,
               log_lik = fit_log_lik,
               AICc = fit_AICc)

}

# Function to combine stats for each protein
combine_model_fit_stats <- function(fit_details){
  agg_fit_details <-
   aggregate(cbind(RSS, n_coeffs, n_obs) ~ Protein_ID,
             data = fit_details,
             FUN = sum)

  agg_fit_details$sigma <- sqrt(agg_fit_details$RSS/agg_fit_details$n_obs)
  agg_fit_details$log_lik <-
    -agg_fit_details$n_obs / 2 * log(2 * pi * agg_fit_details$sigma^2) -
    agg_fit_details$RSS/(2 * agg_fit_details$sigma^2)

  agg_fit_details
}
