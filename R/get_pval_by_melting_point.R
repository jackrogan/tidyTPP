#' Get p-value by comparison of modelled melting point differences
#'
#' @description
#'  Calculate p-value from melting-point differences as described by Savitsky
#'  \emph{et al.} 2014:
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
#' @param fit_tbl A data frame (or `tibble`) of parameters from fitting melting
#'  curves to each protein and taking melting point differences. This must
#'  include the columns:
#'  * `Protein_ID` - Unique protein identifier
#'  * `Condition` - Experimental conditio \emph{e.g.} "Control", "Treatment"
#'  * `min_R_sq` - Minimum coefficient of determination (\eqn{R^2}) for all
#'    observations of a protein
#'  * `max_control_plateau` - Maximum plateau parameter (\eqn{pl}) across all
#'    control curves for a protein
#'  * `diff_melt_point` - The difference in melting point (\eqn{\Delta T_m})
#'    between two observations of a protein, \emph{e.g.} between a treated
#'    and control sample
#'  * `min_slope` - the slope (\eqn{f'(T_{infl})}) of the steepest curve per
#'    comparison
#'
#' @return A `tibble` containing all data in fit_tbl, with additional calculated
#'  column for FDR-adjusted p-value (`adj_pvalue`)
#'
#' @references
#'  Savitski M. M. \emph{et al.} Tracking cancer drugs in living cells by
#'  thermal profiling of the proteome. \emph{Science}, 346: 1255784 (2014)
#'
#'  Benjamini, Y., and Hochberg, Y. Controlling the false discovery
#'  rate: a practical and powerful approach to multiple testing. \emph{Journal
#'  of the Royal Statistical Society Series B}, 57, 289-300 (1995)
#'
#' @export
#'
#' @examples
#' # Minimal data - 4-protein melting curve statistics
#' x <- analysis_data_4prot
#'
#' # Melting point p-value generation
#' get_pval_by_melting_point(x)
#'
get_pval_by_melting_point <- function(fit_tbl){
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

  cols_to_keep <-
    names(pval_tbl) %in% c("Protein_ID", "Condition", "Replicate", "adj_pvalue")
  pval_tbl <-
    pval_tbl[cols_to_keep]

  # Merge tables
  fit_tbl <- merge(fit_tbl, pval_tbl, all.x = TRUE)
  fit_tbl
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
