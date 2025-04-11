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
  # 1. Apply filter criteria - must have 2+ peptide matches
  jointP_tbl <- TPP_tbl[TPP_tbl$Pep_N > 1,]
  # 2. Find jointP - proteins present in all conditions
  jointP_tbl_list <-
    split(jointP_tbl, jointP_tbl[c("Condition", "Replicate")])
  jointP <-
    Reduce(function(x, y) intersect(x, y$Protein_ID),
           jointP_tbl_list[2:length(jointP_tbl_list)],
           init = jointP_tbl_list[[1]]$Protein_ID)
  jointP_tbl_list <-
    lapply(jointP_tbl_list, function(x) x[x$Protein_ID %in% jointP,])

  if(!silent) cat("jointP contains", length(jointP), "Proteins.\n")

  # 3. Filter to fold-change selection
  if(is.null(norm_subset_filter)) {
    norm_subset_filter <- get_normP_filter_tbl(unique(TPP_tbl$Temp))
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

  if(to_plot) {
    median_tbl$Protein_ID <- "Median normP"
    median_tbl$Experiment <- paste(median_tbl$Condition, median_tbl$Replicate)
    plot(build_observed_TPP_plot(median_tbl, facets = FALSE) +
      ggplot2::facet_wrap(ggplot2::vars(.data$Experiment)) +
      ggplot2::ggtitle("Median normP Melting Curve"))
  }

  # 6. Fit curves to median normP sets.
  median_tbl_list <- split(median_tbl, median_tbl[c("Condition", "Replicate")])
  # TODO

  print(median_tbl_list)

  if(!silent) cat("--------------------\n")
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
