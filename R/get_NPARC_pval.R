get_NPARC_pval <- function(TPP_tbl,
                           degrees_of_freedom = 4:7,
                           max_cores = 1,
                           comparisons = NULL,
                           control_name = "Control",
                           to_plot = TRUE,
                           to_save = NULL,
                           silent = FALSE,
                           ...){

  # Check and make sure comparisons are given
  if(is.null(comparisons)){
    conds <- unique(TPP_tbl$Condition)
    reps <- unique(TPP_tbl$Replicate)
    comparisons <- create_comparisons_tbl(conds, reps, "Control")
  }

  comparisons <- comparisons[c("Condition_01", "Condition_02")]
  comparisons <- distinct(comparisons)

  # Loop over condition combinations
  NPARC_tbl <- NULL
  for(k in 1:nrow(comparisons)){
    condition_subset <-
      c(comparisons[[k, "Condition_01"]], comparisons[[k, "Condition_02"]])
    if(!silent) {
      cat("NPARC F-score calculation - ", condition_subset[1], " vs ",
          condition_subset[2], ":\n", sep = "")
    }
    NPARC_tbl_k <- TPP_tbl[TPP_tbl$Condition %in% condition_subset,]

    if(is.numeric(degrees_of_freedom)){
      # Get models per degree of freedom:
      if(!silent) cat("Fitting alternative and null models\n")
      time_total <- proc.time()

      # Parallel:
      if(max_cores > 1){
        n_cores <-
          min(parallel::detectCores(), max_cores, length(degrees_of_freedom))
        cat("Cores:", n_cores, "\nEstimated process time: ",
            round(
              0.21*n_cores +
                (nrow(TPP_tbl)*length(degrees_of_freedom)*0.000113) / n_cores,
              2),
            "s\n")

        cl <- parallel::makeCluster(n_cores)
        parallel::clusterExport(cl, list(
          "get_spline_model_details",
          "get_fit_details",
          "get_natural_spline_lm"
        ),
        envir = environment())

        all_df_splines <-
          parallel::parLapply(cl,
                              1:length(degrees_of_freedom),
                              fit_hypothesis_splines,
                              x_tbl = TPP_tbl,
                              d_f_total = degrees_of_freedom) |>
          Reduce(rbind, x = _)

        parallel::stopCluster(cl)

        # Serially:
      } else {
        # Time taken should be 0.00043 s / Protein row
        cat("Estimated process time: ",
            round(nrow(TPP_tbl) * length(degrees_of_freedom) * 0.000108 , 2),
            "s\n")

        all_df_splines <-
          lapply(1:length(degrees_of_freedom),
                 fit_hypothesis_splines,
                 x = TPP_tbl,
                 d_f_total = degrees_of_freedom) |>
          Reduce(rbind, x = _)
      }
      time_out <- proc.time() - time_total
      if(!silent) cat("Elapsed time:", time_out[["elapsed"]], "s\n")

      if(any(all_df_splines$success)){

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

        # Reshape wide for comparisons
        NPARC_tbl_k <-
          stats::reshape(best_model_splines,
                         direction = "wide",
                         idvar = "Protein_ID",
                         v.names = c("degrees_of_freedom", "RSS", "sigma",
                                     "n_coeffs", "n_obs", "log_lik", "AICc",
                                     "success"),
                         timevar = "hypothesis")

        # Ordinary F-score - as Storey 2005:
        NPARC_tbl_k$F_score <-
          (NPARC_tbl_k$`RSS.H-null` - NPARC_tbl_k$`RSS.H-alt`) /
          NPARC_tbl_k$`RSS.H-alt`

        # Moderate/rescale F-score and p-value:
        # Residual degrees of freedom
        NPARC_tbl_k$`resid_df.H-alt` <-
          NPARC_tbl_k$`n_obs.H-alt` - NPARC_tbl_k$`n_coeffs.H-alt`
        # Posterior variances - using limma squeezeVar
        squeezed_vars <- limma::squeezeVar(NPARC_tbl_k$`sigma.H-alt`^2,
                                           NPARC_tbl_k$`resid_df.H-alt`)
        NPARC_tbl_k$`post_var.H-alt` <- squeezed_vars$var.post
        # Prior distribution degrees of freedom
        NPARC_tbl_k$`prior_df.H-alt` <- squeezed_vars$df.prior
        # Moderated F-score
        NPARC_tbl_k$F_score_mod <-
          (NPARC_tbl_k$`RSS.H-null` - NPARC_tbl_k$`RSS.H-alt`) /
          NPARC_tbl_k$`post_var.H-alt`
        # Degrees of freedom for scaling
        NPARC_tbl_k$num_df <-
          NPARC_tbl_k$`n_coeffs.H-alt` - NPARC_tbl_k$`n_coeffs.H-null`
        NPARC_tbl_k$denom_df <-
          NPARC_tbl_k$`n_obs.H-alt` - NPARC_tbl_k$`n_coeffs.H-alt`
        NPARC_tbl_k$denom_df_adj <-
          NPARC_tbl_k$denom_df + NPARC_tbl_k$`prior_df.H-alt`
        # Scaled F-score
        NPARC_tbl_k$F_scaled <- NPARC_tbl_k$F_score_mod / NPARC_tbl_k$num_df
        # P-value
        NPARC_tbl_k$p_NPARC <-
          stats::pf(NPARC_tbl_k$F_scaled,
                    df1 = NPARC_tbl_k$num_df,
                    df2 = NPARC_tbl_k$denom_df_adj,
                    lower.tail = FALSE)

        # Adjust p-value.
        NPARC_tbl_k$p_adj_NPARC <-
          stats::p.adjust(NPARC_tbl_k$p_NPARC, method = "fdr")

        if(!is.null(to_save) | to_plot){
          title_k <- paste(condition_subset, collapse = " vs ")
          f_plot_k <- build_f_density_plot(NPARC_tbl_k, title_k)
          p_plot_k <- build_p_hist_plot(NPARC_tbl_k, title_k)

          if(to_plot) {
            suppressWarnings(plot(f_plot_k))
            suppressWarnings(plot(p_plot_k))
          }
          if(!is.null(to_save)) {
            save_k <-
              sub("(.*)(\\.[^\\.]+)$",
                  paste0("\\1_", condition_subset[1],
                         "_", condition_subset[2], "_F_score\\2"),
                  to_save)
            suppressWarnings(ggplot2::ggsave(save_k, f_plot_k))
            save_k <-
              sub("(.*)(\\.[^\\.]+)$",
                  paste0("\\1_", condition_subset[1],
                         "_", condition_subset[2], "_p_value\\2"),
                  to_save)
            suppressWarnings(ggplot2::ggsave(save_k, p_plot_k))
          }
        }

        # Keep necessary columns
        NPARC_tbl_k$Condition  <-  condition_subset[1]
        NPARC_tbl_k <-
          NPARC_tbl_k[c("Protein_ID", "Condition", "F_scaled", "p_adj_NPARC")]

      } else {
        if(!silent) {
          cat("No curves successfully fitted:",
              "Try more appropriate degrees of freedom.\n")
        }
      }
    } else {
      if(!silent) {
        cat("Cannot calculate NPARC scores: degrees of freedom not numeric.\n")
      }
    }

    NPARC_tbl <- c(NPARC_tbl, NPARC_tbl_k)
  }

  tibble::as_tibble(NPARC_tbl)
}

# Function to plot F-score density distribution
build_f_density_plot <- function(f_tbl, f_title = NULL){
  f_tbl$facet_label <-
    paste0("df1 = ", f_tbl$num_df,
           "; df2 = ", f_tbl$denom_df)

  f_score_plot <-
    ggplot2::ggplot(f_tbl, ggplot2::aes(x = F_scaled)) +
    ggplot2::geom_density(
      ggplot2::aes(fill = "Calculated"),
      show.legend = FALSE) +
    ggplot2::geom_line(
      ggplot2::aes(
        y = stats::df(F_scaled, df1 = num_df, df2 = denom_df)),
      linewidth = 0.75) +
    ggplot2::theme_bw() +
    ggplot2::lims(x = c(0,10)) +
    ggplot2::facet_wrap(~ facet_label) +
    ggplot2::labs(
      title = f_title,
      subtitle = "F-Score distribution by degrees of freedom",
      x = "Scaled F-statistic",
      y = "Density"
    )
}

# Function to plot p-value density histogram
build_p_hist_plot <- function(p_tbl, p_title = NULL){
  p_val_plot <-
    ggplot2::ggplot(p_tbl, ggplot2::aes(x = p_adj_NPARC)) +
    ggplot2::geom_histogram(
      ggplot2::aes(y = ..density.., fill = "Calculated"),
      boundary = 0,
      bins = 30,
      show.legend = FALSE,
      colour = "black") +
    ggplot2::geom_line(
      ggplot2::aes(y = stats::dunif(p_adj_NPARC)),
      linewidth = 0.75) +
    ggplot2::theme_bw() +
    # ggplot2::facet_wrap(~ facet_label) +
    ggplot2::labs(
      title = p_title,
      subtitle = "p-value distribution",
      x = "Adjusted p-value",
      y = "Density"
    )
}

