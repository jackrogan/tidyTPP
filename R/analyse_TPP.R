#' Analyse TPP-TR Data
#'
#' @description
#' `analyse_TPP()` transforms \emph{Thermal Protein Profiling} (TPP) relative
#'  intensity data, fitting sigmoidal melting curves to the presumably
#'  normalised data, generating curve parameters - including melting points,
#'  then calculating statistical parameters - optionally including melting-point
#'  difference-based p-value as described by Savitsky \emph{et al.} 2014, and
#'  NPARC F-score and p-value as described by Childs \emph{et al.} 2019.
#
#' @inheritParams normalise_TPP
#' @param comparisons A data.frame containing the information needed to build
#'  the comparisons to use to calculate \eqn{\Delta T_m} and other statistics.
#'  Columns must be \emph{Condition_01}, \emph{Replicate_01},
#'  \emph{Condition_02} and \emph{Replicate_01}, \emph{e.g}:
#'
#'  | Condition_01 | Replicate_01 | Condition_02 | Replicate_01 |
#'  | ------------:| ------------:| ------------:| ------------:|
#'  | Treated      | 01           | Control      | 01           |
#'  | Treated      | 02           | Control      | 02           |
#'
#'  If no table is given, one will be generated based on the conditions in the
#'  data, comparing each condition to control and like-for-like replicate
#'  comparisons.
#'
#' @param control_name Character. Character string that matches the control
#'  experiment in `TPP_tbl$Condition`
#' @param p_value_methods Character. Vector of names of statistical significance
#'  approaches to use.
#'
#'  Include "melting_point" to calculate p-value from melting-point differences
#'  as described by Savitsky \emph{et al.} 2014:
#'  * Melting point differences (\eqn{\Delta T_m}) are found and filtered to
#'    those for which \eqn{R^2 > 0.8} for all observations and
#'    \eqn{plateau < 0.3} for all control curves.
#'  * Proteins are ordered by steepness of melting curve slope - steepest first.
#'  * These proteins are then divided into bins of 300, and the final group
#'    added to the second-to-last if less than 300.
#'  * Per bin, the left- and right-sided robust standard deviation is estimated
#'    using the 15.87, 50, and 84.13 percentiles and calculating p-values for
#'    all measurements binwise.
#'  * The p-values thus calculated are adjusted by applying the
#'    Benjamini-Hochberg procedure to control the false discovery rate
#'    (\emph{FDR}).
#'
#'  Include "NPARC" to calculate F-Score and p-value with nonparametric
#'  analysis of response curves \emph{(NPARC)} as described by Childs
#'  \emph{et al.} 2019:
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
#' @param to_plot Boolean. If true, will generate distribution plots of
#'  F-statistic and p-value on calculating \emph{NPARC} statistics
#' @param to_save Character. If supplied, will generate distribution plots of
#'  F-statistic and p-value and save with [ggsave()]
#' @param ... Further arguments to be passed to [fit_melt_by_experiment()] and
#' [nls_multstart()]
#'
#' @return A `tibble` containing all data in TPP_tbl, with additional calculated
#'  columns (where possible) detailing curve parameters and comparison
#'  statistics
#'
#' @references
#'  Savitski M. M. \emph{et al.} Tracking cancer drugs in living cells by
#'  thermal profiling of the proteome. \emph{Science}, 346: 1255784 (2014)
#'
#'  Childs, D., \emph{et al.} Non-Parametric Analysis of Thermal Proteome
#'  Profiles Reveals Novel Drug-Binding Proteins. \emph{Molecular & Cellular
#'  Proteomics}, 18, 2506-2515, (2019)
#'
#' @export
#'
#' @examples
#' # Minimal data - four-protein melt curve
#' x <- quan_data_4prot
#' norm_x <- normalise_TPP(x)
#'
#' # Analyse - single core curve fitting
#' analyse_TPP(norm_x, max_cores = 1)
#'
#' # Analyse - melting-point p-value only
#' analyse_TPP(norm_x, max_cores = 1, p_value_methods = "melting_point")
#'
#' # Analyse - NPARC p-value only and plot F-Score, p-value distributions
#' analyse_TPP(norm_x, max_cores = 1, p_value_methods = "NPARC", to_plot = TRUE)
#'
#' # Custom comparisons, e.g. compare treated replicates
#' comparison_tbl <-  data.frame(
#'   Condition_01 = c("Treated", "Treated", "Treated"),
#'   Replicate_01 = c("01", "02", "01"),
#'   Condition_02 = c("Control", "Control", "Treated") ,
#'   Replicate_02 = c("01", "02", "02")
#' )
#'
#' analyse_TPP(
#'   norm_x,
#'   comparisons = comparison_tbl,
#'   max_cores = 1
#' )
analyse_TPP <-
  function(TPP_tbl,
           comparisons = NULL,
           quantity_column = "rel_quantity",
           control_name = "Control",
           p_value_methods = c("melting_point", "NPARC"),
           silent = FALSE,
           to_plot = FALSE,
           to_save = NULL,
           ...){

  TPP_tbl <- mask_column(TPP_tbl, quantity_column, "quantity")
  if(!"Replicate" %in% colnames(TPP_tbl)) TPP_tbl$Replicate <- "01"

  if(!silent){
    cat("--------------------\n")
    cat("TPP Analysis\n")
    cat("--------------------\n")
    cat("Melting curve fit:\n")
  }

  # Get curve fits
  fit_tbl <-
    fit_melt_by_experiment(TPP_tbl, y_column = "quantity", silent = silent, ...)

  if(!silent){
    cat("--------------------\n")
    cat("Statistic calculation:\n")
  }

  # Default comparisons - assumes condition, replicate columns:
  # non-control vs control, match replicates - TODO  Condition mask
  reps <- unique(TPP_tbl$Replicate)
  if(is.null(comparisons)){
    conds <- unique(TPP_tbl$Condition)
    comparisons <- create_comparisons_tbl(conds, reps, control_name)
  }

  if(!silent){
    cat("\nComparisons:\n")
    print(comparisons)
  }

  # Add control-vs-control
  if(length(reps) > 1){
    controls <- create_control_comparison_tbl(reps, control_name)
    controlled_comparisons <- rbind(comparisons, controls)

    if(!silent){
      cat("\nControl-vs-control:\n")
      print(controls)
    }
  } else {
    if(!silent){
      cat("\nNo control-vs-control stats: not enough replicates\n")
    }
  }

  # Combine with fit data, Get melting Point differences
  if(!silent) cat("\nCalculate per-curve statistics...\n")
  fit_tbl <-
    split(fit_tbl, fit_tbl$Protein_ID) |>
    lapply(find_melting_point_diffs, controlled_comparisons) |>
    Reduce(rbind, x = _)

  # Get statistics per protein: min R2, max vehicle plateau, min slope
  fit_tbl <-
    find_exp_stats(fit_tbl, "min", c("R_sq", "slope"), "Protein_ID")
  control_max_tbl <-
    fit_tbl[fit_tbl$Condition == control_name,] |>
    find_exp_stats("max", "plateau", "Protein_ID")
  control_max_tbl <- control_max_tbl[c("Protein_ID", "max_plateau")]
  control_max_tbl <- unique(control_max_tbl)
  colnames(control_max_tbl)[2] <- "max_control_plateau"
  fit_tbl <- merge(fit_tbl, control_max_tbl)

  if(!silent) cat("Calculate p-values...\n")
  # Get P values (from Tm as in Savitsky 2014)
  if("melting_point" %in% p_value_methods) {
    fit_tbl <- get_pval_by_melting_point(fit_tbl, comparisons)
  }

  # Get NPARC p-values (as in Childs 2019)
  if("NPARC" %in% p_value_methods){
    NPARC_tbl <-
      parse_pass_dots(get_NPARC_pval,
                      add_args = list(TPP_tbl = TPP_tbl,
                                      comparisons = comparisons,
                                      control_name = control_name,
                                      to_plot = to_plot,
                                      to_save = to_save,
                                      silent = silent),
                      ...)

    fit_tbl <- merge(NPARC_tbl, fit_tbl, all = TRUE)
  }

  # Merge tibbles for ease of passing through workflow
  TPP_tbl <- mask_column(TPP_tbl, "quantity", quantity_column)
  full_tbl <- merge(TPP_tbl, fit_tbl,
                    by = c("Protein_ID", "Condition", "Replicate"),
                    all.x = TRUE)

  if(!silent) cat("\nAnalysed.\n--------------------\n")

  tibble::as_tibble(full_tbl)

  }

# Function to create comparison data.frame - all conditions compared to control
create_comparisons_tbl <- function(conds, reps, control_name){
  conds_no_ctl <- conds[conds != control_name]
  comparisons <-
    data.frame(Condition_01 = rep(conds_no_ctl, times = length(reps)),
               Replicate_01 = rep(reps, each = length(conds_no_ctl)),
               Condition_02 =
                 rep(control_name, length(reps) * length(conds_no_ctl)),
               Replicate_02 = rep(reps, each = length(conds_no_ctl)))
}

# Function to create comparison data.frame for control-vs-control stats
create_control_comparison_tbl <- function(reps, control_name){
  control_matrix <- utils::combn(reps, 2)
  controls <-
    data.frame(Condition_01 = control_name,
               Replicate_01 = control_matrix[1,],
               Condition_02 = control_name,
               Replicate_02 = control_matrix[2,])
}

# Function to find melting point differences and min comparison slope
# using comparison data.frame
find_melting_point_diffs <- function(fit_tbl, comparisons){
  comparison_tbl <-
    stats::reshape(comparisons, direction = "long", varying = c(1:4),
                   sep = "_", idvar = "comparison", timevar = "order") |>
    merge(fit_tbl[,c("Protein_ID", "Condition", "Replicate",
                     "melt_point", "slope")]) |>
    stats::reshape(direction = "wide", idvar = "comparison", timevar = "order",
                   v.names = c("Condition", "Replicate", "melt_point", "slope"))

  if("Condition.2" %in% colnames(comparison_tbl)){
    comparison_tbl$comparison <-
      paste(comparison_tbl$Condition.1, comparison_tbl$Replicate.1, "vs",
            comparison_tbl$Condition.2, comparison_tbl$Replicate.2, sep = "_")
    comparison_tbl$diff_melt_point <-
      comparison_tbl$melt_point.1 - comparison_tbl$melt_point.2
    comparison_tbl$min_comparison_slope <-
      pmin(comparison_tbl$slope.1, comparison_tbl$slope.2)
    comparison_tbl <-
      comparison_tbl[,c("Protein_ID", "Condition.1", "Replicate.1", "comparison",
                        "diff_melt_point", "min_comparison_slope")]
    colnames(comparison_tbl) <-
      c("Protein_ID", "Condition", "Replicate", "Comparison",
        "diff_melt_point", "min_comparison_slope")
    comparison_tbl <- merge(fit_tbl, comparison_tbl, all.x = TRUE)
  } else {
    comparison_tbl <- fit_tbl
  }
  tibble::as_tibble(comparison_tbl)
}


