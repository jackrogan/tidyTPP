#' Analyse TPP-TR Data
#'
#' @description
#' `analyse_TPP()` transforms \emph{Thermal Protein Profiling} (TPP) relative
#'  intensity data, fitting sigmoidal melting curves to the presumably
#'  normalised data, generating curve parameters - including melting points,
#'  then calculating statistical parameters including p-value as described by
#'  Savitsky \emph{et al.} 2014:
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
#' @param control_name Character: Character string that matches the control
#'  experiment in `TPP_tbl$Condition`
#' @param ... Further arguments to be passed to [fit_melt_by_experiment()] and
#' [nls_multstart()]
#'
#' @return A `tibble` containing all data in TPP_tbl, with additional calculated
#'  columns (where possible) detaling curve parameters and comparison statistics
#'
#' @references
#'  Savitski M. M. \emph{et al.}, Tracking cancer drugs in living cells by
#'  thermal profiling of the proteome. \emph{Science}, 346: 1255784 (2014)
#'
#'  Benjamini, Y., and Hochberg, Y. Controlling the false discovery
#'  rate: a practical and powerful approach to multiple testing. \emph{Journal
#'  of the Royal Statistical Society Series B}, 57, 289–300 (1995)
#'
#' @export
#'
#' @examples
#' # Minimal data - two-protein melt curve
#' x <- quan_data_normP
#' norm_x <- normalise_TPP(x)
#'
#' # Analyse - single core curve fitting
#' analyse_TPP(norm_x, max_cores = 1)
#'
#' # Custom comparisons, e.g. compare treated replicates
#' comparison_tbl <-  data.frame(
#'   Condition_01 = c("Treated", "Treated", "Treated"),
#'   Replicate_01 = c("01", "02", "01"),
#'   Condition_02 = c("Control", "Control","Treated") ,
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
           silent = FALSE,
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
    comparisons <- rbind(comparisons, controls)

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
    lapply(find_melting_point_diffs, comparisons) |>
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

  # Get p-values (from Tm as in Savitsky 2014)
  # 1. Filter to min R2 > 0.8, max vehicle plateau < 0.3
  pval_tbl <-
    fit_tbl[fit_tbl$min_R_sq > 0.8 & fit_tbl$max_control_plateau < 0.3,]
  pval_tbl <- pval_tbl[!is.na(pval_tbl$diff_melt_point),]

  # 2. Order proteins by ascending min slope (of protein)
  pval_tbl <- pval_tbl[order(pval_tbl$min_slope),]

  # 3. Divide proteins into bins of 300 (plus end bin of 300+)
  initial_bins <- (nrow(pval_tbl) %/% 300) + 1
  pval_tbl$bin <-
    rep(1:initial_bins, each = 300, length.out = nrow(pval_tbl))
  if(nrow(pval_tbl[pval_tbl$bin == initial_bins,]) < 300){
    pval_tbl[pval_tbl$bin == initial_bins, "bin"] <- initial_bins - 1
  }
  pval_tbl_list <- split(pval_tbl, pval_tbl$bin)

  # 4. Per bin, estimate the left- and right-sided robust standard deviation
  #    using the 15.87, 50, and 84.13 percentiles and calculating
  #    P values for all measurements
  pval_tbl_list <- lapply(pval_tbl_list, calculate_binwise_pvalue)
  pval_tbl <- Reduce(rbind, pval_tbl_list)

  # 5. Adjust with benjamini-hochberg over full sample set.
  pval_tbl$adj_pvalue <- stats::p.adjust(pval_tbl$pvalue, "BH")

  pval_tbl <- pval_tbl[c("Protein_ID", "Condition", "Replicate",
                         "Comparison", "adj_pvalue")]

  # Merge tibbles for ease of pipeline and export
  fit_tbl <- merge(fit_tbl, pval_tbl, all.x = TRUE)
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

# Function to find melting point differences with comparison data.frame
find_melting_point_diffs <- function(fit_tbl, comparisons){
  comparison_tbl <-
    stats::reshape(comparisons, direction = "long", varying = c(1:4),
            sep = "_", idvar = "comparison", timevar = "order") |>
    merge(fit_tbl[,c("Protein_ID", "Condition", "Replicate",
                     "melt_point")]) |>
    stats::reshape(direction = "wide", idvar = "comparison", timevar = "order",
            v.names = c("Condition", "Replicate", "melt_point"))

  comparison_tbl$comparison <-
    paste(comparison_tbl$Condition.1, comparison_tbl$Replicate.1, "vs",
          comparison_tbl$Condition.2, comparison_tbl$Replicate.2, sep = "_")
  comparison_tbl$diff_melt_point <-
    comparison_tbl$melt_point.2 - comparison_tbl$melt_point.1

  comparison_tbl <-
    comparison_tbl[,c("Protein_ID", "Condition.1", "Replicate.1", "comparison",
                      "diff_melt_point")]
  colnames(comparison_tbl) <-
    c("Protein_ID", "Condition", "Replicate", "Comparison", "diff_melt_point")
  comparison_tbl <- merge(fit_tbl, comparison_tbl, all.x = TRUE)

  tibble::as_tibble(comparison_tbl)
}

# Function to compute per-bin p-values
calculate_binwise_pvalue <- function(binned_data){
  # Use 15.87, 50, and 84.13 percentiles
  mp_quantiles <-
    stats::quantile(binned_data$diff_melt_point,
             probs = c(0.1587, 0.5, 0.8413),
             na.rm=TRUE)
  q1 <- mp_quantiles[1]
  q2 <- mp_quantiles[2]
  q3 <- mp_quantiles[3]

  # Separate left- and right-sided data
  right_binned_data <- binned_data[binned_data$diff_melt_point > q2,]
  left_binned_data <- binned_data[binned_data$diff_melt_point <= q2,]

  # Generate z-scores
  right_binned_data$z <-  (right_binned_data$diff_melt_point - q2) / (q3-q2)
  left_binned_data$z <- (q2 - left_binned_data$diff_melt_point) / (q2-q1)

  # Recombine, calculate p-value
  binned_data <- rbind(right_binned_data, left_binned_data)
  binned_data$pvalue <- 2 * stats::pnorm(abs(binned_data$z), lower.tail = FALSE)
  binned_data <- binned_data[order(binned_data$min_slope),]

}
