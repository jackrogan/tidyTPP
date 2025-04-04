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
#'
#' @param data TPP protein observation data; see [build_melt_plot_with_model()]
#' @param fit_parameters parameters to predict protein melting curve; see
#'  [predict_melt_curve()]
#' @param to_plot Boolean. If true, run [plot()] on `ggplot` object
#' @param to_save Character. If supplied, save plot with [ggsave()]
#' @param to_add_to_ggplot List. Optional list of further additions to `ggplot`;
#'  will be added using [+.gg()]
#' @param ... Additional arguments to be parsed to [predict_melt_curve()] or
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
#' # Customise plot
#' plot_melt(
#'   x,
#'   x_parameters,
#'   to_add_to_ggplot = ggtitle("ATIC", subtitle = "Protein Melting Curve"),
#'   n_predict = 200,
#'   annotate = "both",
#'   facets = FALSE)
#'
plot_melt <- function(data,
                      fit_parameters,
                      to_plot = TRUE,
                      to_save = NULL,
                      to_add_to_ggplot = NULL,
                      ...
                      ){

  # Deal with arguments to pass to functions
  parse_pass_dots <- function(func, add_args = list(), ...){
    dots <- list(...)
    func_dots <- dots[names(dots) %in% names(formals(func))]
    do.call(func, c(add_args, func_dots))
  }

  model_data <-
    parse_pass_dots(predict_melt_curve,
                    list("fit_parameters" = fit_parameters),
                    ...)
  melting_plot <-
    parse_pass_dots(build_melt_plot_with_model,
                    list("data" = data, "predicted_data" = model_data),
                    ...)

  if(!is.null(to_add_to_ggplot)){
    addition_list <- c(list(melting_plot), list(to_add_to_ggplot))
    melting_plot <- Reduce("+", addition_list)
  }
  if(to_plot) plot(melting_plot)
  if(!is.null(to_save)) {
    ggsave(to_save, melting_plot)
  }

  data
}
