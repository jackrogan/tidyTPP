# Normalise TPP data according to Savitsky 2014
# Return (normalised) tidy tibble

normalise_TPP <- function(TPP_tbl,
                          to_plot = FALSE,
                          to_save = NULL,
                          quality_filter = NULL,
                          norm_subset_filter = NULL,
                          quantity_column = NULL,
                          silent = FALSE){

  if(!silent){
    cat("--------------------\n")
    cat("TPP Normalisation\n")
    cat("--------------------\n")
  }

  if(is.null(quality_filter)) quality_filter <- get_npeptide_quality_filter()
  print(quality_filter)
  TPP_temps <- unique(TPP_tbl$Temp)
  # 1. Apply filter criteria - must have 2+ peptide matches
  TPP_tbl <- TPP_tbl[TPP_tbl$Pep_N > 1,]
  # 2. Find jointP - proteins present in all conditions
  jointP_tbl_list <-
    split(TPP_tbl, TPP_tbl[c("Condition", "Replicate")])
  jointP <-
    Reduce(function(x, y) intersect(x, y$Protein_ID),
           jointP_tbl_list[2:length(jointP_tbl_list)],
           init = jointP_tbl_list[[1]]$Protein_ID)
  jointP_tbl_list <-
    lapply(jointP_tbl_list, function(x) x[x$Protein_ID %in% jointP,])

  if(!silent) cat("jointP contains", length(jointP), "Proteins.\n")

  # 3. Filter to fold-change selection
  if(is.null(norm_subset_filter)) {
    norm_subset_filter <- get_normP_filter_tbl(TPP_temps)
  }

  if(!silent) {
    cat("\nnormP criteria:\n")
    print(norm_subset_filter)
  }

  normP_list <-
    lapply(jointP_tbl_list, filter_to_normP, norm_subset_filter)

  # 4. Take subset from condition/replicate pair with longest protein list,
  # and apply to all experiments
  which_subset <-
    sapply(normP_list, length) |>
    which.max()
  normP <- normP_list[[which_subset]]

  if(!silent) cat("\nnormP contains", length(normP), "Proteins.\n")

  # 5. Get Median fold-changes for each Temp, Condition, Replicate normP subset
  median_tbl <-
    lapply(jointP_tbl_list, function(x) x[x$Protein_ID %in% normP,]) |>
    Reduce(rbind, x = _) |>
    stats::aggregate(rel_quantity ~ Condition + Replicate + Temp,
                     FUN = stats::median)

  # 6. Fit curves to median normP sets.
  if(!silent) {
    cat("--------------------\n")
    cat("Fit melting curve to normP medians:\n")
  }
  norm_fit_tbl <-
    fit_melt_by_experiment(median_tbl,
                           experiment_cols = c("Condition", "Replicate"),
                           max_cores = 1)

  # Plot normalisation curves
  if(to_plot | !is.null(to_save)) {
    median_tbl$Protein_ID <- "Median normP"
    median_tbl$Exp <-
      paste(median_tbl$Condition, median_tbl$Replicate)
    norm_fit_model_tbl <- predict_melt_curve(norm_fit_tbl)
    norm_fit_model_tbl$Exp <-
      paste(norm_fit_model_tbl$Condition, norm_fit_model_tbl$Replicate)
    norm_melt_plot <-
      build_melt_plot_with_model(median_tbl,
                                 norm_fit_model_tbl,
                                 facets_column = "Exp",
                                 annotate = "R_sq",
                                 rules = FALSE) +
      ggplot2::ggtitle("Median normP Melting Curve")

    if(to_plot) plot(norm_melt_plot)
    if(!is.null(to_save)) ggplot2::ggsave(to_save, norm_melt_plot)
  }

  # 7. Use condition with best fitted curve (by R2)
  best_norm_fit <-
    norm_fit_tbl[which.max(norm_fit_tbl$R_sq),]
  if(!silent){
    cat("\nBest fitted normP median curve:\n")
    norm_fit_out <-
      as.data.frame(best_norm_fit[,c("Condition", "Replicate", "R_sq")])
    norm_fit_out$R_sq <- round(norm_fit_out$R_sq, 3)
    print(norm_fit_out)
  }

  # 8. Find normalisation coefficients
  model_median_tbl <- predict_melt_curve(best_norm_fit, T_seq = TPP_temps)
  model_median_tbl <- model_median_tbl[,c("Temp", "model_quantity")]
  median_tbl <- merge(median_tbl, model_median_tbl)
  median_tbl$norm_coefficient <-
    median_tbl$model_quantity / median_tbl$rel_quantity
  coeff_tbl <- median_tbl[,c("Condition", "Replicate", "Temp", "norm_coefficient")]


  # 9. Normalise all proteins
  TPP_tbl <- merge(TPP_tbl, coeff_tbl)
  TPP_tbl$rel_quantity <- TPP_tbl$rel_quantity * TPP_tbl$norm_coefficient

  if(!silent){
    cat("--------------------\n")
    cat(nrow(TPP_tbl) / nrow(coeff_tbl), "proteins normalised\n")
    cat("--------------------\n")
  }

  tibble::as_tibble(TPP_tbl)
}

# Find closest temperatures used to standard normP temperature filter
get_normP_filter_tbl <- function(data_temps = NULL){
  if(!is.null(data_temps)){
    used_temps <-
    c(data_temps[which.min(abs(data_temps - 56))],
      data_temps[which.min(abs(data_temps - 63))],
      data_temps[which.min(abs(data_temps - 67))])
  } else {
    used_temps <- c(56, 63, 67)
  }

  filter_tbl <-
    data.frame(Temp = used_temps,
               lower = c(0.4, -Inf, -Inf),
               higher = c (0.6, 0.3, 0.2))
}

get_npeptide_quality_filter <- function(){
  data.frame(col = "Pep_N", lower = 2, upper = Inf)
}

# Reduce data to list of proteins that fulfil all criteria
filter_to_normP <- function(x_tbl, filter){
  filter_subset_proteins <- function(sub_no){
    x_tbl[(x_tbl$Temp == filter$Temp[sub_no] &
             x_tbl$rel_quantity > filter$lower[sub_no] &
             x_tbl$rel_quantity < filter$higher[sub_no]),] |>
      _$Protein_ID |>
      unique()
  }
  x_normP <-
    lapply(1:nrow(filter), filter_subset_proteins) |>
    Reduce(intersect, x = _)
}
