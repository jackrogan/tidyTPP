#' Get p-value by nonparametric analysis of response curves \emph{(NPARC)}
#'
#' @description
#'  Calculate F-Score and p-value with nonparametric analysis of response curves
#'  \emph{(NPARC)} as described by Childs \emph{et al.} 2019:
#'  * For each protein and condition comparison, assume a null hypothesis
#'    \eqn{H_{null}} (measurements from both conditions can be modelled with the
#'    same sigmoidal melting curve) and alternate hypothesis \eqn{H_{alt}}
#'    (measurements each condition can be modelled with separate sigmoidal
#'    melting curve.)
#'  * Curves are fitted to each hypothesis for each protein, and the degrees
#'    of freedom estimated for each curve fit.
#'  * The F-statistic and p-value is computed using the residual sum of squares
#'    \emph{RSS} from \eqn{H_{alt}} and \eqn{H_{null}} calculated degrees
#'    of freedom for for each protein and condition comparison.
#'  * p-values are adjusted by applying the
#'    Benjamini-Hochberg procedure to control the false discovery rate
#'    (\emph{FDR}).
#'
#' @inheritParams analyse_TPP
#' @param degrees_of_freedom Numeric. Range of degrees of freedom to use to
#'  generate possible \eqn{H_{null}} and \eqn{H_{alt}} curves.
#' @param all_stats Boolean. If true, report RSS for \eqn{H_{alt}},
#'  \eqn{H_{null}}, and degrees of freedom with output.
#' @param alt_model A `data.frame` with model parameter data for the alternative
#'  hypothesis. Must include experimental condition and replicate columns as
#'  well as `a`, `b`, `plateau`, `melt_point`, `infl_point`, `slope`, `R_sq`,
#'  `RSS`, `sigma`, `n_ceoffs`, `n_obs`, `log_lik`, and `AICc`.
#'
#'  See example output for [fit_melt_by_experiment]
#' @param max_cores Integer. The maximum number of cores to parallelise the
#'  spline fitting over, if less than the maximum number of cores available.
#'  If `max_cores` is 1, then spline fitting will be run serially in
#'  a single process
#'
#' @return A `tibble`, with the columns `Protein_ID`, `Condition`, `F-scaled`,
#'  and `p_adj_NPARC`
#'
#' @references
#'  Childs, D., \emph{et al.} Non-Parametric Analysis of Thermal Proteome
#'  Profiles Reveals Novel Drug-Binding Proteins. \emph{Molecular & Cellular
#'  Proteomics}, 18, 2506-2515, (2019)
#'
#'  Benjamini, Y., and Hochberg, Y. Controlling the false discovery
#'  rate: a practical and powerful approach to multiple testing. \emph{Journal
#'  of the Royal Statistical Society Series B}, 57, 289-300 (1995)
#'
#' @export
#'
#' @examples
#' # Minimal data - two-protein melt curve
#' x <- quan_data_20prot
#' norm_x <- normalise_TPP(x)
#'
#' # NPARC F-score and p-value generation
#' get_NPARC_pval(norm_x,
#'                quantity_column = "rel_quantity",
#'                max_cores = 1)
#'
#' # Plot F-score and p-value distribution
#' get_NPARC_pval(norm_x,
#'                quantity_column = "rel_quantity",
#'                max_cores = 1,
#'                to_plot = TRUE)
get_NPARC_pval <- function(TPP_tbl,
                           NPARC_fit_method = "nls",
                           degrees_of_freedom = c(4:7),
                           alt_model = NULL,
                           max_cores = 8,
                           comparisons = NULL,
                           quantity_column = "quantity",
                           control_name = "Control",
                           all_stats = FALSE,
                           to_plot = FALSE,
                           to_save = NULL,
                           silent = FALSE){

  TPP_tbl <- mask_column(TPP_tbl, quantity_column, "quantity")

  # Check that quantity is given
  if(!"quantity" %in% colnames(TPP_tbl)){
    NPARC_tbl <- NULL

  } else {
    # Check and make sure comparisons are given
    if(is.null(comparisons)){
      conds <- unique(TPP_tbl$Condition)
      reps <- unique(TPP_tbl$Replicate)
      comparisons <- create_comparisons_tbl(conds, reps, control_name)
    }

    comparisons <- comparisons[c("Condition_01", "Condition_02")]
    comparisons <- unique(comparisons)
    print(comparisons)

    # Loop over condition combinations
    NPARC_tbl <- NULL

    for(k in 1:nrow(comparisons)){
      condition_subset <-
        c(comparisons[[k, "Condition_01"]], comparisons[[k, "Condition_02"]])
      if(!silent) {
        cat("\nNPARC F-score calculation - ", condition_subset[1], " vs ",
            condition_subset[2], ":\n", sep = "")
      }
      NPARC_tbl_k <- TPP_tbl[TPP_tbl$Condition %in% condition_subset,]

      if(NPARC_fit_method == "splines"){
        model_details_tbl <-
          fit_splines_by_df(NPARC_tbl_k,
                            degrees_of_freedom,
                            max_cores = max_cores,
                            silent = silent)
      } else {
        if(!silent & !NPARC_fit_method %in% c("nonlinear", "nls")) {
          cat("Fit method not recognised - using non-linear model\n")
        }

        NPARC_tbl_k <- mask_column(NPARC_tbl_k, "quantity", "rel_quantity")

        # Use already-calculated sigmoid curve fits for alternative hypothesis
        if(!is.null(alt_model)) {
          alt_details_tbl <- alt_model
        } else {
          alt_details_tbl <-
            fit_melt_by_experiment(NPARC_tbl_k,
                                   experiment_cols = c("Protein_ID", "Condition"),
                                   silent = silent,
                                   max_cores = max_cores)
        }

        alt_details_tbl <- combine_model_fit_stats(alt_details_tbl)

        alt_details_tbl$hypothesis <- "H-alt"
        alt_details_tbl$AICc <- NA

        alt_details_tbl <-
          alt_details_tbl[c("Protein_ID", "RSS", "sigma", "n_coeffs", "n_obs",
                            "log_lik", "AICc", "hypothesis")]

        null_details_tbl <-
          fit_melt_by_experiment(NPARC_tbl_k,
                                 experiment_cols = c("Protein_ID"),
                                 silent = silent,
                                 max_cores = max_cores)

        null_details_tbl$hypothesis <- "H-null"

        null_details_tbl <-
          null_details_tbl[c("Protein_ID", "RSS", "sigma", "n_coeffs", "n_obs",
                             "log_lik", "AICc", "hypothesis")]

        model_details_tbl <-
          as.data.frame(rbind(alt_details_tbl, null_details_tbl))
        model_details_tbl$degrees_of_freedom <- NA
        model_details_tbl$success <- TRUE
      }

      if(!is.null(model_details_tbl)){

        # Reshape wide for comparisons
        NPARC_tbl_k <-
          stats::reshape(model_details_tbl,
                         direction = "wide",
                         idvar = "Protein_ID",
                         v.names = c("degrees_of_freedom", "RSS", "sigma",
                                     "n_coeffs", "n_obs", "log_lik", "AICc",
                                     "success"),
                         timevar = "hypothesis")

        get_f_score_TPP <- function(NPARC_tbl_k){

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
              suppressWarnings(suppressMessages(plot(f_plot_k)))
              suppressWarnings(suppressMessages(plot(p_plot_k)))
            }
            if(!is.null(to_save)) {
              save_k <-
                sub("(.*)(\\.[^\\.]+)$",
                    paste0("\\1_", condition_subset[1],
                           "_", condition_subset[2], "_F_score\\2"),
                    to_save)
              suppressWarnings(
                suppressMessages(ggplot2::ggsave(save_k, f_plot_k)))
              save_k <-
                sub("(.*)(\\.[^\\.]+)$",
                    paste0("\\1_", condition_subset[1],
                           "_", condition_subset[2], "_p_value\\2"),
                    to_save)
              suppressWarnings(suppressMessages(
                ggplot2::ggsave(save_k, p_plot_k)))
            }
          }

          # Keep necessary columns
          NPARC_tbl_k$Condition  <-  condition_subset[1]
          NPARC_tbl_k <-
            NPARC_tbl_k[c("Protein_ID", "Condition", "num_df", "denom_df",
                          "F_scaled", "p_adj_NPARC")]
          colnames(NPARC_tbl_k) <- c("Protein_ID", "Condition", "d1", "d2",
                                     "F_scaled", "p_adj_NPARC")
          NPARC_tbl_k
        }

        get_f_score_NPARC <- function(NPARC_tbl_k){

          # RSS Difference
          NPARC_tbl_k$RSS.diff <- NPARC_tbl_k$`RSS.H-null` - NPARC_tbl_k$`RSS.H-alt`
          NPARC_tbl_k <- NPARC_tbl_k[!is.na(NPARC_tbl_k$RSS.diff),]
          NPARC_tbl_k <- NPARC_tbl_k[NPARC_tbl_k$RSS.diff > 0,]
          # Estimate degrees of freedom empirically:
          # Estimate sigma0 squared with median and median absolute deviation
          NPARC_tbl_k$sigma0_sq <-
            0.5 *
            stats::mad(NPARC_tbl_k$RSS.diff, na.rm = TRUE)^2 /
            stats::median(NPARC_tbl_k$RSS.diff, na.rm = TRUE)

          # Fit Chi2 estimate for degrees of freedom
          NPARC_tbl_k$num_df <-
            MASS::fitdistr(
              x =  NPARC_tbl_k$RSS.diff/ NPARC_tbl_k$sigma0_sq,
              densfun = "chi-squared",
              start = list(df = 1),
              method = "Brent",
              lower = 0,
              upper = nrow(NPARC_tbl_k)
            )[["estimate"]]

          NPARC_tbl_k$denom_df <-
            MASS::fitdistr(
              x =  NPARC_tbl_k$`RSS.H-alt`/ NPARC_tbl_k$sigma0_sq,
              densfun = "chi-squared",
              start = list(df = 1),
              method = "Brent",
              lower = 0,
              upper = nrow(NPARC_tbl_k)
            )[["estimate"]]

          # NPARC_tbl_k$`RSS.H-null` <- NPARC_tbl_k$`RSS.H-null`/ NPARC_tbl_k$sigma0_sq
          # NPARC_tbl_k$`RSS.H-alt` <- NPARC_tbl_k$`RSS.H-alt`/ NPARC_tbl_k$sigma0_sq
          # NPARC_tbl_k$`RSS.diff` <- NPARC_tbl_k$`RSS.diff`/ NPARC_tbl_k$sigma0_sq

          # Scaled F-statistic
          NPARC_tbl_k$F_scaled <-
            (NPARC_tbl_k$RSS.diff / (NPARC_tbl_k$num_df * NPARC_tbl_k$sigma0_sq)) /
            (NPARC_tbl_k$`RSS.H-alt` / (NPARC_tbl_k$denom_df * NPARC_tbl_k$sigma0_sq))

          # Adjusted p-value
          NPARC_tbl_k$pvalue <-
            1 - stats::pf(NPARC_tbl_k$F_scaled,
                          df1 = NPARC_tbl_k$num_df,
                          df2 = NPARC_tbl_k$denom_df)
          NPARC_tbl_k$p_adj_NPARC <-  stats::p.adjust(NPARC_tbl_k$pvalue, "BH")

          if(!is.null(to_save) | to_plot){
            title_k <- paste(condition_subset, collapse = " vs ")
            f_plot_k <- build_f_density_plot(NPARC_tbl_k, title_k)
            p_plot_k <- build_p_hist_plot(NPARC_tbl_k, title_k)

            if(to_plot) {
              suppressWarnings(suppressMessages(plot(f_plot_k)))
              suppressWarnings(suppressMessages(plot(p_plot_k)))
            }
            if(!is.null(to_save)) {
              save_k <-
                sub("(.*)(\\.[^\\.]+)$",
                    paste0("\\1_", condition_subset[1],
                           "_", condition_subset[2], "_F_score\\2"),
                    to_save)
              suppressWarnings(
                suppressMessages(ggplot2::ggsave(save_k, f_plot_k)))
              save_k <-
                sub("(.*)(\\.[^\\.]+)$",
                    paste0("\\1_", condition_subset[1],
                           "_", condition_subset[2], "_p_value\\2"),
                    to_save)
              suppressWarnings(suppressMessages(
                ggplot2::ggsave(save_k, p_plot_k)))
            }
          }
          # Keep necessary columns
          NPARC_tbl_k$Condition  <-  condition_subset[1]
          if(all_stats){
            NPARC_tbl_k <-
              NPARC_tbl_k[c("Protein_ID", "Condition", "RSS.diff", "RSS.H-alt",
                            "RSS.H-null", "num_df", "denom_df",
                            "F_scaled", "p_adj_NPARC")]
            colnames(NPARC_tbl_k) <- c("Protein_ID", "Condition", "RSS_diff",
                                       "RSS_H_alt", "RSS_H_null", "d1", "d2",
                                       "F_scaled", "p_adj_NPARC")
          } else {
            NPARC_tbl_k <-
              NPARC_tbl_k[c("Protein_ID", "Condition",
                            "F_scaled", "p_adj_NPARC")]
            colnames(NPARC_tbl_k) <- c("Protein_ID", "Condition",
                                       "F_scaled", "p_adj_NPARC")
          }

          NPARC_tbl_k
        }

        NPARC_tbl_k <- get_f_score_NPARC(NPARC_tbl_k)
        NPARC_tbl <- rbind(NPARC_tbl, NPARC_tbl_k)
      }
    }
  }

  if(!is.null(NPARC_tbl)){
    NPARC_tbl <- mask_column(NPARC_tbl, "quantity", quantity_column)
    tibble::as_tibble(NPARC_tbl)
  } else {
    TPP_tbl <- mask_column(TPP_tbl, "quantity", quantity_column)
    TPP_tbl
  }
}

# Plot F-score density distribution
#' @importFrom rlang .data
build_f_density_plot <- function(f_tbl, f_title = NULL){
  f_tbl$facet_label <-
    paste0("df1 = ", signif(f_tbl$num_df, 3),
           "; df2 = ", signif(f_tbl$denom_df, 3))

  f_score_plot <-
    ggplot2::ggplot(f_tbl, ggplot2::aes(x = .data$F_scaled)) +
    ggplot2::geom_density(
      ggplot2::aes(fill = "Calculated"),
      show.legend = FALSE) +
    ggplot2::geom_line(
      ggplot2::aes(
        y = stats::df(.data$F_scaled,
                      df1 = .data$num_df,
                      df2 = .data$denom_df)),
      linewidth = 0.75) +
    ggplot2::theme_bw() +
    ggplot2::lims(x = c(0,10)) +
    ggplot2::facet_wrap(~ .data$facet_label) +
    ggplot2::labs(
      title = f_title,
      subtitle = "F-Score distribution by degrees of freedom",
      x = "Scaled F-statistic",
      y = "Density"
    )
}

# Plot p-value density histogram
#' @importFrom rlang .data
build_p_hist_plot <- function(p_tbl, p_title = NULL){
  p_val_plot <-
    ggplot2::ggplot(p_tbl, ggplot2::aes(x = .data$p_adj_NPARC)) +
    ggplot2::geom_histogram(
      ggplot2::aes(
        y = ggplot2::after_stat(!!str2lang("density")), fill = "Calculated"),
      boundary = 0,
      bins = 30,
      show.legend = FALSE,
      colour = "black") +
    ggplot2::geom_line(
      ggplot2::aes(y = stats::dunif(.data$p_adj_NPARC)),
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

