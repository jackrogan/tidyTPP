# Get hits (or top n hits), including call to plot

get_TPP_hits <- function(
    TPP_data,
    n_hits = NULL,
    hit_criteria = list(
      pvalue_min_threshold = 0.05,
      # pvalue_min_threshold = 0.25,
      pvalue_max_threshold = 0.1,
      # pvalue_max_threshold = 1.1,
      DTm_same_sign = TRUE,
      DTm_gt_Dcontrol = TRUE,
      slope_threshold = -0.06
    ),
    control_name = "Control",
    to_export = "TPP_hits.xlsx",
    to_plot = TRUE,
    to_save = NULL,
    plot_separately = FALSE,
    silent = FALSE,
    ...
){
  TPP_hits <- TPP_data

  if(!silent){
    cat("--------------------\n")
    cat("TPP Hit Identification\n")
    cat("--------------------\n")

    cat("Hit Criteria:\n")
    print.data.frame(
      data.frame(hit_criteria),
      row.names = FALSE)
    cat("--------------------\n")
  }

  TPP_hits <-
    TPP_hits[, colnames(TPP_hits) %in%
               c("Protein_ID", "Condition", "Replicate", "Comparison", "a", "b",
                 "plateau", "melt_point", "infl_point", "slope", "R_sq",
                 "diff_melt_point", "min_comparison_slope",
                 "min_R_sq", "min_slope", "max_control_plateau",
                 "adj_pvalue")]
  TPP_hits <- unique(TPP_hits)

  # Move control delta-Tm to its own column.
  control_tbl <- TPP_hits[TPP_hits$Condition == control_name,]

  # Use mean absolute control difference
  control_tbl <-
    aggregate(cbind(melt_point, diff_melt_point) ~ Protein_ID,
              control_tbl,
              FUN = \(x) mean(abs(x), na.rm = TRUE),
              na.action = na.pass)
  colnames(control_tbl)[2] <- "mean_control_melt_point"
  colnames(control_tbl)[3] <- "abs_diff_melt_control"
  TPP_hits <- merge(TPP_hits[TPP_hits$Condition != control_name,], control_tbl)

  starting_colnames <- colnames(TPP_hits)
  stat_columns_present <- c("adj_pvalue", "diff_melt_point",
                            "min_comparison_slope")
  stat_columns_present <-
    stat_columns_present[which(stat_columns_present %in% colnames(TPP_hits))]

  TPP_hits <-
    find_exp_stats(TPP_hits,
                    stat_func = c("min", "max", "all_same_sign", "min_abs"),
                    stat_column = stat_columns_present,
                    experiment_cols = c("Protein_ID", "Condition"))

  if(length(hit_criteria) > 0){

    # Filter 1: pvalue threshold
    if("pvalue_min_threshold" %in% names(hit_criteria)){
      TPP_hits <-
          TPP_hits[TPP_hits$min_adj_pvalue < hit_criteria$pvalue_min_threshold,]

    }
    if("pvalue_max_threshold" %in% names(hit_criteria)){
        TPP_hits <-
          TPP_hits[TPP_hits$max_adj_pvalue < hit_criteria$pvalue_max_threshold,]
    }

    # Filter 2: Both melting point differences have the same sign
    if("DTm_same_sign" %in% names(hit_criteria)) {
      TPP_hits <- TPP_hits[TPP_hits$all_same_sign_diff_melt_point,]
    }

    # Filter 3: Minimum absolute melting point difference greater than
    # mean absolute control difference
    if("DTm_gt_Dcontrol" %in% names(hit_criteria)) {
      TPP_hits <-
        TPP_hits[
          TPP_hits$min_abs_diff_melt_point > TPP_hits$abs_diff_melt_control,]
    }

    # Filter 4: steeper slope in comparison pair below threshold
    # Allow not filtering if comparison minimum slope column not present
    if("slope_threshold" %in% names(hit_criteria)){
      if("max_min_comparison_slope" %in% colnames(TPP_hits)){
        TPP_hits <-
          TPP_hits[
            TPP_hits$max_min_comparison_slope < hit_criteria$slope_threshold,]
      } else {
        if(!silent) cat("Warning: Cannot filter by steepest comparison slope,",
                        "missing 'min_comparison_slope' column.\n")
      }
    }
    TPP_hits <- TPP_hits[c(starting_colnames, "max_adj_pvalue")]
  }

  # Summary table of hits
  TPP_hits$Comparison <-
    gsub("^(.*)_\\d+(_vs_.*)_\\d+$", "\\1\\2", TPP_hits$Comparison)

  TPP_hits <-
    aggregate(
      cbind(max_adj_pvalue, melt_point, mean_control_melt_point,
            diff_melt_point, abs_diff_melt_control) ~
        Protein_ID + Condition + Comparison,
      TPP_hits,
      FUN = mean)
  colnames(TPP_hits) <-
    c("Protein_ID", "Condition", "Comparison", "max_adj_pvalue", "mean_melt_point",
      "mean_control_melt_point", "mean_diff_melt_point",
      "mean_control_diff_melt_point")


  # Get top n hits in order of ascending pvalue
  TPP_hits <- TPP_hits[order(TPP_hits$max_adj_pvalue),]
  if(!is.null(n_hits)) TPP_hits <- utils::head(TPP_hits, n_hits)

  # Export data
  if(!is.null(to_export)) {
    if(!silent) cat("Exporting hit data...\n")
    parse_pass_dots(export_TPP,
                    list("TPP_data" = TPP_hits, "file_name" = to_export),
                    ...)
  }

  # Plot top hits
  if(to_plot | !is.null(to_save)){
    if(!silent) cat("Plotting hit melting curves...\n")
    hit_proteins <- unique(TPP_hits$Protein_ID)
    hit_quan_data <- TPP_data[which(TPP_data$Protein_ID %in% hit_proteins),]
    if(plot_separately){
      hit_quan_data <- split(hit_quan_data, ~ Protein_ID)
    } else {
      hit_quan_data <- list(hit_quan_data)
    }
    lapply(hit_quan_data, plot_melt, to_plot = to_plot, to_save = to_save, ...)
  }

  tibble::as_tibble(TPP_hits)
  if(!silent) cat("--------------------\n")
}
