get_pval_by_melting_point <- function(fit_tbl){
  # Get p-values (from Tm as in Savitsky 2014) - using means of replicates
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

  pval_tbl <- pval_tbl[c("Protein_ID", "Condition", "Replicate", "adj_pvalue")]

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
