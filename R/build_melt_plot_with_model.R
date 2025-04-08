#' Build Protein Melt Plot Using Predicted Data
#'
#' @description
#' `build_melt_plot_with_model()` creates a [ggplot2] plot, showing the protein
#'  melting curves given in `data` and `predicted_data` for observed points and
#'  fitted line, respectively.
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
#'
#' @param predicted_data A data frame (or [tibble]) containing predicted protein
#'  melting curve data \emph{e.g.} from [predict_melt_curve()]. As for `data`,
#'  this must contain the columns \emph{Temp} and \emph{rel_quantity} (or name
#'  specified with `quan_name`), and the columns \emph{Protein_ID},
#'  \emph{Condition} and \emph{Replicate} will be used for `ggplot` aesthetics,
#'  but defaults will be supplied if they are missing.
#' @param annotate Control which text annotations should be added to the plot:
#'  * "R_sq" adds \eqn{R^2} fit regression values (top-right)
#'  * "melt_point" adds \eqn{T_m} melting points (bottom-left)
#'  * "both" and "none" add both and none, respectively
#' @param rules Boolean. Whether to add rules showing \eqn{T_m}
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
#' x_predicted <- predict_melt_curve(MP_data_ATIC)
#'
#' # Create plot using default appearance
#' build_melt_plot_with_model(x, x_predicted)
#'
#' # Specify column names
#' build_melt_plot_with_model(
#'   data = x,
#'   predicted_data = x_predicted,
#'   quan_column = "rel_quantity",
#'   colour_column = "Condition",
#'   shape_column = "Replicate"
#'   )
#'
#' # Remove annotations and facets
#' build_melt_plot_with_model(
#'   data = x,
#'   predicted_data = x_predicted,
#'   annotate = "none",
#'   rules = FALSE,
#'   facets = FALSE
#'   )
#'
#' # As a ggplot object, additions can be made to returned plot
#' # For example, the default behaviour is to hide the "shape" legend,
#' # however, this can be disabled:
#' build_melt_plot_with_model(x, x_predicted) +
#'   ggplot2::scale_shape_discrete()

build_melt_plot_with_model <-
  function(data,
           predicted_data,
           annotate = c("R_sq", "melt_point", "both", "none"),
           rules = TRUE,
           facets = TRUE,
           quan_column = "rel_quantity",
           colour_column = "Condition",
           shape_column = "Replicate"
           ){

    # Select default annotation behaviour
    if(any(sum(annotate == c("R_sq", "melt_point", "both", "none")) == c(0,4))){
      annotate <- "none"
    }

    # Check for lack of condition/replicate columns
    check_experiment_cols <- function(tbl){
      if(!any(colnames(tbl) == "Protein_ID")) tbl$Protein_ID <- "Protein"
      if(!any(colnames(tbl) == "Condition")) tbl$Condition <- tbl$Protein_ID
      if(!any(colnames(tbl) == "Replicate")) tbl$Replicate <- "01"
      tbl
    }

    data <- check_experiment_cols(data)
    predicted_data <- check_experiment_cols(predicted_data)

    # Set quantity columns
    if(!any(colnames(data) == quan_column)) {
      quan_column <- "rel_quantity"
    }
    colnames(predicted_data)[
      colnames(predicted_data) == "model_quantity"] <- "quantity"

  # Plot observed data, add predicted curves
  melt_plot <-
    build_observed_TPP_plot(data, facets,
                            quan_column,
                            colour_column,
                            shape_column) +
    ggplot2::geom_line(data = predicted_data)

  # Add rules on plot for melting points if required (default yes)
  if(rules){
    melt_plot <-
      melt_plot +
      ggplot2::geom_segment(data = predicted_data,
                            ggplot2::aes(x = min(.data$Temp),
                                         xend = .data$melt_point,
                                         y = 0.5,
                                         yend = 0.5),
                            linetype = "dashed", colour = "grey") +
      ggplot2::geom_segment(data = predicted_data,
                            ggplot2::aes(x = .data$melt_point,
                                         xend = .data$melt_point,
                                         y = 0,
                                         yend = 0.5),
                            linetype = "dashed")
  }

  # Add chosen annotations (default none)
  if(any(annotate == c("R_sq", "melt_point", "both"))){
    max_q <- max(data$quantity)
    annotation_data <-
      predicted_data[c("Condition", "Replicate", "R_sq", "melt_point")]
    annotation_data <- unique(annotation_data)
    annotation_data$Experiment <-
      as.factor(paste0(annotation_data$Condition,
                       '"~"',
                       annotation_data$Replicate))
    annotation_data$Exp_no <- as.integer(annotation_data$Experiment)

    if(any(annotate == c("R_sq", "both"))){
      annotation_data$R_sq_label <-
        paste0('R^2~("', annotation_data$Experiment,
               '"):~', sprintf('"%1.3f"', signif(annotation_data$R_sq, 3)))
      annotation_data$R_sq_y <- max_q * (1 - (annotation_data$Exp_no - 1) / 15)

      melt_plot <-
        melt_plot +
        ggplot2::geom_text(data = annotation_data,
                           ggplot2::aes(x = max(data$Temp),
                                        y = .data$R_sq_y,
                                        label = .data$R_sq_label),
                           parse = TRUE,
                           hjust = 1,
                           show.legend = FALSE)
    }

    if(any(annotate == c("melt_point", "both"))){
      annotation_data$Tm_label <-
        paste0('T[m]~("', annotation_data$Experiment,
               '"):~', sprintf('"%2.1f"', signif(annotation_data$melt_point, 3)),
               "~degree*C")
      annotation_data$Tm_y <-
        max_q * (max(annotation_data$Exp_no) - annotation_data$Exp_no) / 15

      melt_plot <-
        melt_plot +
        ggplot2::geom_text(data = annotation_data,
                           ggplot2::aes(x = min(data$Temp),
                                        y = .data$Tm_y,
                                        label = .data$Tm_label),
                           parse = TRUE,
                           hjust = 0,
                           vjust = 0,
                           show.legend = FALSE)
    }
  }

  melt_plot
}
