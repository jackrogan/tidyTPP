# Function to loop over all degrees of freedom and fit splines
fit_splines_by_df <-
  function(TPP_subtbl, degrees_of_freedom, max_cores, silent = FALSE){
  if(is.numeric(degrees_of_freedom)){
    # Get models per degree of freedom:
    if(!silent) cat("Fitting alternative and null models\n")
    time_total <- proc.time()

    # Parallel:
    if(max_cores > 1){
      n_cores <-
        min(parallel::detectCores(), max_cores, length(degrees_of_freedom))
      if(!silent){
        est_core_time <- nrow(TPP_subtbl) * length(degrees_of_freedom) * 0.000113
        est_time <- 0.21 * n_cores + est_core_time / n_cores
        cat("Cores:", n_cores,
            "\nEstimated process time: ", round(est_time, 2), "s\n")
      }
      cl <- parallel::makeCluster(n_cores)
      parallel::clusterExport(cl, list(
        "get_model_fit_stats",
        "get_fit_details",
        "get_natural_spline_lm"
      ),
      envir = environment())

      all_df_splines <-
        parallel::parLapply(cl,
                            1:length(degrees_of_freedom),
                            fit_hypothesis_splines,
                            x_tbl = TPP_subtbl,
                            d_f_total = degrees_of_freedom,
                            silent = silent) |>
        Reduce(rbind, x = _)

      parallel::stopCluster(cl)

      # Serially:
    } else {
      if(!silent){
        # Time taken should be 0.00043 s / Protein row
        est_time <- nrow(TPP_subtbl) * length(degrees_of_freedom) * 0.000108
        cat("Estimated process time: ", round(est_time , 2), "s\n")
      }
      all_df_splines <-
        lapply(1:length(degrees_of_freedom),
               fit_hypothesis_splines,
               x = TPP_subtbl,
               d_f_total = degrees_of_freedom,
               silent = silent) |>
        Reduce(rbind, x = _)
    }
    time_out <- proc.time() - time_total
    if(!silent) cat("Elapsed time:", time_out[["elapsed"]], "s\n")

    if(any(all_df_splines$success & all_df_splines$hypothesis == "H-alt")){

      # Select degrees of freedom that minimise AICc for each protein
      # break ties with minimum complexity
      best_model_splines <- all_df_splines[all_df_splines$success,]
      best_model_splines <-
        best_model_splines[all_df_splines$hypothesis == "H-alt",]
      best_model_splines <-
        best_model_splines[order(best_model_splines$Protein_ID,
                                 best_model_splines$AICc,
                                 best_model_splines$degrees_of_freedom),]
      best_model_splines <-
        best_model_splines[!duplicated(best_model_splines$Protein_ID),]

      # filter all spline tbl to best models
      best_model_splines <-
        best_model_splines[c('Protein_ID', 'degrees_of_freedom')]
      best_model_splines <-
        merge(best_model_splines, all_df_splines)

    } else {
      if(!silent) {
        cat("No alternate hypothesis curves successfully fitted:",
            "Try more appropriate degrees of freedom.\n")
        best_model_splines <- NULL
      }
    }
  } else {
    if(!silent) {
      cat("Cannot calculate NPARC scores: degrees of freedom not numeric.\n")
      best_model_splines <- NULL
    }
  }

  best_model_splines
}

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
  lm_details <- get_model_fit_stats(spline_lm)

  cbind(data.frame(Protein_ID = x_tbl$Protein_ID[1]), lm_details)
}
