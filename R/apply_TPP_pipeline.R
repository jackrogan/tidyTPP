#' Apply tidyTPP TPP-TR pipeline with default settings
#'
#' @description
#' Apply tidyTPP TPP-TR pipeline with default settings: by default imports
#' Spectronaut report \emph{.csv} file.
#'
#' @inheritParams import_custom_format
#' @param to_plot Boolean. Whether to generate plots to the console:
#'  * Normalisation curves
#'  * NPARC analysis F-score and p-value distribution
#'  * Hit protein melting curves
#' @param max_cores Integer. Maximum number of cores to run analysis using. If
#'  1, will use no parallel or multi-core functions.
#'
#' @return A `tibble` of hit proteins, sorted by NPARC p-value
#' @export
#'
#' @examples
#'  # Minimal Example - ATIC and CFL1 data
#'  normP_report_file <-
#'    system.file("extdata", "normP_report.csv", package = "tidyTPP")
#'  config_file <-
#'    system.file("extdata", "ATIC_config.csv", package = "tidyTPP")
#'
#'  # Analyse using single core, generate and show plots
#'  normP_hits <-
#'    apply_TPP_pipeline(datafile = normP_report_file,
#'                       config = config_file,
#'                       to_plot = TRUE,
#'                       max_cores = 1)
#'
#'  normP_hits
apply_TPP_pipeline <- function(datafile,
                               config,
                               path = NULL,
                               to_plot = FALSE,
                               max_cores = 4){

  # (include "/" if not given)
  if(!is.null(path)){
    path <- sub("([^/])$", "\\1/", path)
    full_path <- paste0(path, datafile)
  } else {
    full_path <- datafile
  }

  import_spectronaut(datafile, config, path) |>

  normalise_TPP(to_plot = to_plot) |>

  analyse_TPP(to_plot = to_plot, max_cores = max_cores) |>

  export_TPP(sub("\\.[^\\.]+$", "_Results.xlsx", full_path)) |>

  get_TPP_hits(to_export = sub("\\.[^\\.]+$", "_Hits.xlsx", full_path),
               to_plot = to_plot)

}
