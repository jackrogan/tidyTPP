#' Apply tidyTPP TPP-TR pipeline with default settings
#'
#' @description
#' Apply tidyTPP TPP-TR pipeline with default settings: by default imports
#' Spectronaut report \emph{.csv} file.
#'
#' @inheritParams import_TPP
#' @param import_format Character. Specify a format to import with an
#'  \emph{import_} function:
#'
#'  * "Spectronaut" or "SN" - peptide report from \emph{Spectronaut} DIA protein
#'  quantification (\emph{.csv}): [import_spectronaut]
#'
#'  * "ProteomeDiscover" or "PD" - exported proteins table from
#'  \emph{Thermo Proteome Discoverer} DDA quantification (\emph{.txt}):
#'  [import_proteomediscoverer]
#'
#'  * For a custom format import, use [import_TPP] directly
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
#'  # Minimal Example - four-protein data
#'  temp_loc <- tempdir()
#'  source_report <-
#'    system.file("extdata", "4_protein_peptide_report.csv", package = "tidyTPP")
#'  source_config <-
#'    system.file("extdata", "4_protein_config.csv", package = "tidyTPP")
#'  file.copy(c(source_report, source_config), temp_loc)
#'
#'  # Analyse using single core, generate and show plots
#'  normP_hits <-
#'    apply_TPP_pipeline(datafile = "4_protein_peptide_report.csv",
#'                       config = "4_protein_config.csv",
#'                       path = temp_loc,
#'                       import_format = "Spectronaut",
#'                       to_plot = TRUE,
#'                       max_cores = 1)
#'
#'  unlink(
#'   c(
#'     paste0(temp_loc, "\\", "4_protein_peptide_report.csv"),
#'     paste0(temp_loc, "\\", "4_protein_config.csv"),
#'     paste0(temp_loc, "\\", "4_protein_peptide_Results.xlsx"),
#'     paste0(temp_loc, "\\", "4_protein_peptide_Hits.xlsx"))
#'   )
#'
#'  normP_hits
apply_TPP_pipeline <- function(datafile,
                               config,
                               path = NULL,
                               import_format = "Spectronaut",
                               to_plot = FALSE,
                               max_cores = 4){

  # (include "/" if not given)
  if(!is.null(path)){
    path <- sub("([^/])$", "\\1/", path)
    full_path <- paste0(path, datafile)
  } else {
    full_path <- datafile
  }

  import_TPP(datafile, config, path, import_format) |>

  normalise_TPP(to_plot = to_plot) |>

  analyse_TPP(to_plot = to_plot, max_cores = max_cores) |>

  export_TPP(sub("\\.[^\\.]+$", "_Results.xlsx", full_path)) |>

  get_TPP_hits(to_export = sub("\\.[^\\.]+$", "_Hits.xlsx", full_path),
               to_plot = to_plot)

}
