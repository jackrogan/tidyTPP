#' Get Hit Protein Results from TPP-TR data
#'
#' @description
#' `get_TPP_hits()` filters TPP-TR result data to candidates for thermally
#'  stabilised or destabilised proteins, and returns a per-protein summary of
#'  results. Default hit criteria areas described by Savitsky et al. 2014:
#'  * The melting point difference (\eqn{\Delta T_m}) compared against control
#'    for a protein has a minimum  \emph{FDR}-adjusted p-value of less than 0.05
#'    and a maximum of less than 0.10 across biological replicates.
#'  * All replicates show a positive \eqn{\Delta T_m} or all negative.
#'  * All \eqn{\Delta T_m} values compared against controls are greater than the
#'    maximum \eqn{\Delta T_m} between controls.
#'  * For each comparison \emph{versus} control, the steepest slope is steeper
#'    than -0.06.
#'
#'  Melting point curves for hits will be generated if `to_plot` is `TRUE`
#'  or `to_save` is not `NULL`.
#'
#' @inheritParams analyse_TPP
#' @param TPP_data TPP-TR Analysis results data from `analyse_TPP`. Must have
#' `melt_point`, `diff_melt_point`, and `adj_pvalue` columns.
#'
#' `min_comparison_slope` (the steepest slope parameter between the plots
#'  compared for each \eqn{\Delta T_m} comparison) is necessary for the
#' `slope_threshold` hit criteria to be in effect.
#' @param n_hits Optional integer. Maximum number of hits to be retained in
#'  summary or exported. Hits are ranked by ascending adjusted p-value, and the
#'  top `n_hits` kept.
#' @param hit_criteria List of criteria for selecting hits. Default criteria:
#'
#'  |                     |        |
#'  | --------------------:| -----:|
#'  | pvalue_min_threshold | 0.05  |
#'  | pvalue_max_threshold | 0.10  |
#'  | DTm_same_sign        | TRUE  |
#'  | DTm_gt_Dcontrol      | TRUE  |
#'  | slope_threshold      | -0.6  |
#'
#'  * \emph{pvalue_min_threshold:} Threshold \emph{at least one} adjusted
#'    p-value must meet per protein
#'  * \emph{pvalue_max_threshold:} Threshold \emph{every} adjusted p-value must
#'    meet per protein
#'  * \emph{DTm_same_sign:} \eqn{\Delta T_m} must have the same sign
#'  * \emph{DTm_gt_Dcontrol:} \eqn{\Delta T_m} must be
#'    \eqn{> \Delta T_m (controls)}
#'  * \emph{slope_threshold:} Steepness threshold for comparison steepest slope.
#'
#' @param to_export Character. File path to save hit data. If `NULL`,
#'  this will not be saved.
#' @param to_plot Boolean. If true, will generate melting point curves of hits
#'  and `plot()`
#' @param to_save Character. If supplied, will generate melting point curves of
#'  hits and save with [ggsave()]
#' @param plot_separately Boolean. If TRUE, and either `to_plot` is TRUE or
#'  `to_save` is not `NULL`, melting point plots will be generated per protein.
#'  If false, this will be a plot with proteins as facets by default.
#' @param ... Additional arguments to be passed to [export_TPP()],
#' [predict_melt_curve()] or [build_melt_plot_with_model()]
#'
#' @return A `tibble` containing columns for adjusted pvalue, treated and
#' control mean \eqn{T_m} and \eqn{\Delta T_m}. One row per protein/condition
#' combination.
#'
#' @references
#'  Savitski M. M. \emph{et al.}, Tracking cancer drugs in living cells by
#'  thermal profiling of the proteome. \emph{Science}, 346: 1255784 (2014)
#'
#' @export
#'
#' @examples
#' # Minimal data - analysed single-protein dataset
#' x <- MP_stat_data_ATIC
#' hit_results_file <- tempfile(pattern = "hits_", fileext = ".xlsx")
#'
#' # Default criteria return no hits - p-values do not meet criteria
#' hits <- get_TPP_hits(x, to_export = NULL, to_plot = FALSE)
#' hits
#'
#' # Using only subset of hit criteria returns single protein to export as .xlsx
#' hits <-
#'   get_TPP_hits(
#'     x,
#'     hit_criteria = list(DTm_same_sign = TRUE,
#'                         DTm_gt_Dcontrol = TRUE),
#'     to_export = hit_results_file,
#'     to_plot = FALSE)
#' hits
#'
#' unlink(hit_results_file)
get_TPP_hits <- function(
    TPP_data,
    n_hits = NULL,
    hit_criteria = list(
      pvalue_min_threshold = 0.05,
      pvalue_max_threshold = 0.1,
      DTm_same_sign = TRUE,
      DTm_gt_Dcontrol = TRUE,
      slope_threshold = -0.06
    ),
    control_name = "Control",
    to_export = "TPP_hits.xlsx",
    to_plot = FALSE,
    to_save = NULL,
    plot_separately = FALSE,
    silent = FALSE,
    ...
){
  TPP_hits <- TPP_data
  if(!"Replicate" %in% colnames(TPP_hits)) TPP_hits$Replicate <- "01"

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
                 "adj_pvalue")]
  TPP_hits <- unique(TPP_hits)

  # Move control delta-Tm to its own column.
  control_tbl <- TPP_hits[TPP_hits$Condition == control_name,]

  # Use mean absolute control difference
  control_tbl <-
    stats::aggregate(cbind(melt_point, diff_melt_point) ~ Protein_ID,
                     control_tbl,
                     FUN = \(x) mean(abs(x), na.rm = TRUE),
                     na.action = stats::na.pass)
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
        if(!silent) {
          message("Warning: Cannot filter by steepest comparison slope, ",
                  "missing 'min_comparison_slope' column.\n")
        }
      }
    }
    TPP_hits <- TPP_hits[c(starting_colnames, "max_adj_pvalue")]
  }

  if(nrow(TPP_hits) > 0 & !all(is.na(TPP_hits$Protein_ID))){

    # Summary table of hits
    TPP_hits$Comparison <-
      gsub("^(.*)_\\d+(_vs_.*)_\\d+$", "\\1\\2", TPP_hits$Comparison)

    TPP_hits <-
      stats::aggregate(
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
    if((to_plot | !is.null(to_save)) & any(colnames(TPP_data) == "Temp")){
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

    if(!silent) {
      cat(length(unique(TPP_hits$Protein_ID)),
          "hits found.\n--------------------\n")
    }
  } else {
    if(!silent) cat("No hits found.\n--------------------\n")
    TPP_hits <- TPP_hits[!is.na(TPP_hits$Protein_ID),]
  }

  tibble::as_tibble(TPP_hits)
}
