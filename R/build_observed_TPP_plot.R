#' Build Plot from TPP observed Protein data
#'
#' @description
#' `build_observed_TPP_plot()` creates a [ggplot2] plot, showing the protein
#'  melting data points given in `data` only. To add fitted melting curves use
#'  `build_melt_plot_with_model()`.
#'
#' @param data A data frame (or [tibble]) containing proteomics data from a
#'  thermal protein profiling (TPP) experiment. In all cases, this must contain
#'  the columns:
#'  \describe{
#'  \item{Temp}{Experiment temperature in \eqn{\degree C}, will be used for the x
#'  axis values}
#'  \item{rel_quantity}{(or name specified with `quan_name`) Protein quantity
#'  measurement corresponding to `Fraction non-denatured`, will be used for the
#'  y-axis values}
#'  }
#'  The following columns will be used for the `ggplot` `aesthetics`, but
#'  default values will be assumed if not present:
#'  \describe{
#'  \item{Protein_ID}{Gene name identifying the protein - used for facets and
#'  to provide values for default `Condition` column; defaults to generic
#'  "Protein"}
#'  \item{Condition}{Treatment given to this sample set - used for colour;
#'  defaults to Protein_ID values}
#'  \item{Replicate}{Replicate number for this sample set - used for shape;
#'  defaults to "01"}
#'  }
#' @param facets Boolean. Whether to separate facets using \emph{Protein_ID}
#' @param quan_column Character. The column in `data` and `predicted_data` used
#' for plot y-values.
#' @param colour_column Character. The column in `data` and `predicted_data` used
#' for plot colour aesthetic.
#' @param shape_column Character. The column in `data` and `predicted_data` used
#' for plot shape aesthetic.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @importFrom rlang .data
#'
#' @examples
#' # Minimal data - ATIC melt curve
#' x <- quan_data_ATIC
#'
#' # Create plot using default appearance
#' build_observed_TPP_plot(x)
#'
#' # Specify column names
#' build_observed_TPP_plot(
#'   data = x,
#'   quan_column = "rel_quantity",
#'   colour_column = "Condition",
#'   shape_column = "Replicate"
#'   )
#'
#' # Remove facets
#' build_observed_TPP_plot(
#'   data = x,
#'   facets = FALSE
#'   )
#'
#' # As a ggplot object, additions can be made to returned plot
#' # For example, changing the theme:
#' build_observed_TPP_plot(x) +
#'   ggplot2::theme_minimal()
build_observed_TPP_plot <- function(data,
                                    facets = TRUE,
                                    quan_column = "rel_quantity",
                                    colour_column = "Condition",
                                    shape_column = "Replicate"
){
  # Plot observed TPP data Points
  colnames(data)[colnames(data) == quan_column] <- "quantity"

  # Deal with missing or present columns
  get_aesthetic_column <- function(aes_name){
    if(any(colnames(data) == aes_name)) {
      aes_name <- data[[aes_name]]
    } else{
      aes_name <- NULL
    }
  }
  colour_column_aes <- get_aesthetic_column(colour_column)
  shape_column_aes <- get_aesthetic_column(shape_column)

  melt_plot <-
    ggplot2::ggplot(data,
                    ggplot2::aes(x = .data$Temp, y = .data$quantity,
                                 colour = colour_column_aes,
                                 shape = shape_column_aes)) +
    ggplot2::geom_point() +
    ggplot2::labs(title = "Protein Melting Curve",
                  x = "Temperature / \u00B0C",
                  y = "Fraction Non-denatured",
                  colour = colour_column,
                  shape = shape_column) +
    ggplot2::scale_shape_discrete(guide = "none")
  if(is.null(colour_column_aes)) {
    melt_plot <- melt_plot + ggplot2::scale_colour_discrete(guide = "none")
  }
  melt_plot <- melt_plot + ggplot2::theme_bw()

  # Facet separate proteins if specified (default yes)
  if(facets){
    melt_plot <-
      melt_plot +
      ggplot2::facet_wrap(ggplot2::vars(.data$Protein_ID))
  }

  melt_plot
}
