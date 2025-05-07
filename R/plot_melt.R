#' Plot Melt Curve
#'
#' @description
#' `plot_melt` builds and plots a protein melting point from proteomics
#'  \emph{Thermal Protein Profiling} (TPP) data and corresponding fitted model.
#'  Generates prediction data with [predict_melt_curve()], and plot with
#'  [build_melt_plot_with_model()] Included options are to plot in an
#'  interactive session, or save using [ggplot2::ggsave()], while returning the
#'  original data set, to match the behaviour of \emph{e.g.} [print()] or
#'  [plot()].
#'
#' @param data TPP protein observation data; see [build_melt_plot_with_model()].
#'  If fit_parameters are not given separately they will be expected as part of
#'  this dataset, as `a`, `b`, `plateau` and, if intending to annotate,
#'  `melt_point`, and `R_sq` columns
#' @param fit_parameters parameters to predict protein melting curve; see
#'  [predict_melt_curve()]. Overrides columns in `data`
#' @param to_plot Boolean. If true, run [plot()] on `ggplot` object
#' @param to_save Character. If supplied, save plot with [ggsave()]
#' @param file_suffix_column Character. The name of a column to be used as an
#'  additional file suffix before the extension if `to_save` is given and the
#'  plot is to be saved. The first element of the given column of `data` will be
#'  used to identify the plot.
#' @param to_add_to_ggplot List. Optional list of further additions to `ggplot`;
#'  will be added using [+.gg()]
#' @param ... Additional arguments to be passed to [predict_melt_curve()] or
#' [build_melt_plot_with_model()]
#'
#' @return A data frame or [tibble] given in `data`
#' @export
#'
#' @examples
#' # Minimal data - ATIC melt curve
#' x <- quan_data_ATIC
#' x_parameters <- MP_data_ATIC
#'
#' # Plot melting point curve
#' plot_melt(x, x_parameters)
#'
#' # Use combined data table
#' x_plus_params <- merge(x, x_parameters)
#' plot_melt(x_plus_params)
#'
#' # Customise plot
#' plot_melt(
#'   x,
#'   x_parameters,
#'   to_add_to_ggplot = ggplot2::ggtitle("ATIC", subtitle = "Protein Melting Curve"),
#'   n_predict = 200,
#'   annotate = "both",
#'   facets = FALSE)
#'
plot_melt <- function(data,
                      fit_parameters = NULL,
                      to_plot = TRUE,
                      to_save = NULL,
                      file_suffix_column = NULL,
                      to_add_to_ggplot = NULL,
                      ...
                      ){

  # Deal with fit parameters - if not given, try to extract from data
  if(is.null(fit_parameters)){
    fit_parameters <-
      data[colnames(data) %in% c("Protein_ID",
                                 "Condition",
                                 "Replicate",
                                 "a", "b", "plateau",
                                 "melt_point",
                                 "R_sq")]
  }

  if(ncol(fit_parameters) > 3){
    model_data <-
      parse_pass_dots(predict_melt_curve,
                      list("fit_parameters" = fit_parameters),
                      ...)
    melting_plot <-
      parse_pass_dots(build_melt_plot_with_model,
                      list("data" = data, "predicted_data" = model_data),
                      ...)
  } else {
    melting_plot <-
      parse_pass_dots(build_observed_TPP_plot,
                      list("data" = data),
                      ...)
  }


  if(!is.null(to_add_to_ggplot)){
    addition_list <- c(list(melting_plot), list(to_add_to_ggplot))
    melting_plot <- Reduce("+", addition_list)
  }

  if(to_plot) plot(melting_plot)
  if(!is.null(to_save)) {
    file_suffix <- NULL
    if(!is.null(file_suffix_column)) {
      file_suffix = paste0("_", data[[1, file_suffix_column]])
    }
    ggplot2::ggsave(
      sub("(\\.[[:alnum:]]+)*$", paste0(file_suffix, "\\1"), to_save),
      melting_plot)
  }

  data
}
