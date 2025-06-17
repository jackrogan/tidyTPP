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
#'  temp_loc <- tempdir()
#'  source_report <-
#'    system.file("extdata", "normP_report.csv", package = "tidyTPP")
#'  source_config <-
#'    system.file("extdata", "ATIC_config.csv", package = "tidyTPP")
#'  normP_report_file <- "normP_report.csv"
#'  config_file <- "ATIC_config.csv"
#'
#'  # Analyse using single core, generate and show plots
#'  normP_hits <-
#'    apply_TPP_pipeline(datafile = normP_report_file,
#'                       config = config_file,
#'                       path = temp_loc,
#'                       to_plot = TRUE,
#'                       max_cores = 1)
#'
#'  unlink(
#'   c(
#'     paste0(temp_loc, "\\", "normP_report.csv"),
#'     paste0(temp_loc, "\\", "ATIC_config.csv"),
#'     paste0(temp_loc, "\\", "normP_report_Results.xlsx"),
#'     paste0(temp_loc, "\\", "normP_report_Hits.xlsx"))
#'   )
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
