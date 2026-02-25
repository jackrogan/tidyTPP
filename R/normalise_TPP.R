#' Normalise TPP-TR Data
#'
#' @description
#'  `normalise_TPP()` transforms \emph{Thermal Protein Profiling} (TPP) relative
#'  intensity data, normalising against fitted median melting curves, as
#'  described by Savitsky \emph{et al.} 2014:
#'  * The set of proteins observed in all experimental groups, \emph{jointP}, is
#'    identified.
#'  * A subset of proteins is identified, \emph{normP}, with a clear \eqn{T_m}
#'    ~ \eqn{56 \degree C} for each experimental group, and the largest
#'    \emph{normP} used.
#'  * Median fold-change values over all proteins in \emph{normP} are calculated
#'    for each experiment, at each temperature point.
#'  * Melting curves are fitted, and the curve with the best \eqn{R^2} used as
#'    the normal model.
#'  * A correction factor is calculated for each median fold-change in order to
#'    coincide with the best-fitting curve.
#'  * These correction factors are applied to all the protein measurements in the
#'    respective experimental groups.
#'
#' @inheritParams plot_melt
#' @inheritParams fit_melting_curve
#' @param TPP_tbl A data frame (or [tibble]) containing proteomics data from a
#'  thermal protein profiling (TPP) experiment. This must contain the columns:
#'  * `Protein_ID`: Unique protein identity
#'  * `Condition` : Category of treated or control sample
#'  * `Replicate` : Replicate number.
#'  * `Temperature`: Experiment temperature in \eqn{\degree C}
#'  * Protein quantity measurement corresponding to the non-denatured fraction
#'
#'  If replicates are not given, a single replicate per condition will be
#'  assumed.
#' @param quality_filter A data.frame containing the information needed to build
#' an initial quantity filter. Columns must be \emph{col}, \emph{lower} and
#' \emph{upper}. Default:
#'
#'  | col          | lower         | upper        |
#'  | ------------ | -------------:| ------------:|
#'  | Pep_N        | 2             | Inf          |
#'
#'  Resulting in filtering for number of detected peptides per protein
#'  \emph{Pep_N} \eqn{\geq 2}
#' @param norm_subset_filter A data.frame containing the information needed to
#'  build the protein filter for normP. Columns must be \emph{Temp}, \emph{lower}
#'  and \emph{upper}, \emph{e.g}:
#'
#'  | Temp         | lower         | upper        |
#'  | ------------:| -------------:| ------------:|
#'  | 56           | 0.4           | 0.6          |
#'  | 63           | -Inf          | 0.3          |
#'  | 67           | -Inf          | 0.2          |
#'
#'  If no table is given, one will be generated based on the one above, using
#'  the closest temperature points to 56, 63, and 67 \eqn{\degree C},
#'  respectively
#' @param quantity_column Character. Name for the column containing
#' protein quantity values in `TPP_tbl`
#'
#' @return A `tibble`, as given in TPP_tbl, with the quantity values
#' (from `quantity_column`) replaced with normalised quantities, and an
#' additional `norm_coefficients` column
#'
#' @references
#'  Savitski M. M. \emph{et al.}, Tracking cancer drugs in living cells by
#'  thermal profiling of the proteome. \emph{Science}, 346: 1255784 (2014)
#'
#' @export
#'
#' @examples
#' # Minimal data - four-protein melt curve
#' x <- quan_data_4prot
#'
#' # Normalise data
#' normalise_TPP(x)
#'
#' # Plot normalisation curves
#' normalise_TPP(x, to_plot = TRUE)
#'
#' # Filter proteins to >= 3 observed peptides
#' normalise_TPP(
#'   x,
#'   quality_filter = data.frame(col = "Pep_N", lower = 3, upper = Inf)
#' )
#'


normalise_TPP <- function(TPP_tbl,
                          to_plot = FALSE,
                          to_save = NULL,
                          quality_filter = NULL,
                          norm_subset_filter = NULL,
                          quantity_column = "rel_quantity",
                          silent = FALSE){
  TPP_temps <- unique(TPP_tbl$Temp)
  if(!silent){
    cat("--------------------\n")
    cat("TPP Normalisation\n")
    cat("--------------------\n")
  }
  if(is.null(quality_filter)) quality_filter <- get_npeptide_quality_filter()
  if(!silent){
    cat("Quality Criteria:\n")
    print(quality_filter)
  }
  TPP_tbl <- mask_column(TPP_tbl, quantity_column, "quantity")
  if(!"Replicate" %in% colnames(TPP_tbl)) TPP_tbl$Replicate <- "01"

  # 1. Apply filter criteria - must have 2+ peptide matches
  if(nrow(quality_filter) > 0){
    for(i in 1:nrow(quality_filter)){
      TPP_tbl <- filter_from_criterion(TPP_tbl, quality_filter[i,])
    }
  }

  # 2. Find jointP - proteins present in all conditions
  jointP_tbl_list <-
    split(TPP_tbl, TPP_tbl[c("Condition", "Replicate")])

  jointP <-
    Reduce(function(x, y) intersect(x, y$Protein_ID),
           jointP_tbl_list[2:length(jointP_tbl_list)],
           init = jointP_tbl_list[[1]]$Protein_ID)
  jointP_tbl_list <-
    lapply(jointP_tbl_list, function(x) x[x$Protein_ID %in% jointP,])

  if(!silent) cat("\njointP contains", length(jointP), "Proteins.\n")

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
    stats::aggregate(quantity ~ Condition + Replicate + Temp,
                     FUN = stats::median)
  # 6. Fit curves to median normP sets.
  if(!silent) {
    cat("--------------------\n")
    cat("Fit melting curve to normP medians:\n")
  }
  norm_fit_tbl <-
    fit_melt_by_experiment(median_tbl,
                           experiment_cols = c("Condition", "Replicate"),
                           y_column = "quantity",
                           max_cores = 1,
                           silent = silent)

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
    median_tbl$model_quantity / median_tbl$quantity
  coeff_tbl <- median_tbl[,c("Condition", "Replicate", "Temp", "norm_coefficient")]


  # 9. Normalise all proteins
  TPP_tbl <- merge(TPP_tbl, coeff_tbl)
  TPP_tbl$quantity <- TPP_tbl$quantity * TPP_tbl$norm_coefficient
  n_prots <- length(unique(TPP_tbl$Protein_ID))

  if(!silent){
    cat("--------------------\n")
    cat(n_prots, "proteins normalised\n")
    cat("--------------------\n")
  }

  TPP_tbl <- mask_column(TPP_tbl, "quantity", quantity_column)
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
               upper = c (0.6, 0.3, 0.2))
}

get_npeptide_quality_filter <- function(){
  data.frame(col = "Pep_N", lower = 2, upper = Inf)
}

# Reduce data to list of proteins that fulfil all criteria
filter_to_normP <- function(x_tbl, filter){
  filter_subset_proteins <- function(sub_no){
    x_tbl[(x_tbl$Temp == filter$Temp[sub_no] &
             x_tbl$quantity > filter$lower[sub_no] &
             x_tbl$quantity < filter$upper[sub_no]),] |>
      _$Protein_ID |>
      unique()
  }
  x_normP <-
    lapply(1:nrow(filter), filter_subset_proteins) |>
    Reduce(intersect, x = _)
}

# Filter by individual quality filter table row
filter_from_criterion <- function(x_tbl, filter_row){
  x_tbl[(x_tbl[,filter_row$col[[1]]] >= filter_row$lower[[1]] &
           x_tbl[,filter_row$col[[1]]] < filter_row$upper[[1]]),]
}
