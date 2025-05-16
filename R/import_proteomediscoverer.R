#' Import TPP data from Proteome Discoverer exports into tidier tibble format.
#'
#' @inheritParams import_custom_format
#' @param datafile Character. Specify a Proteome Discoverer exported .txt
#' filename containing the \emph{Proteins} table.
#'
#' @return A tibble giving the TPP protein quantity data relative to the lowest
#'  temperature, \eqn{T_1}, as \emph{rel_quantity}, in a tidier format: 1 row
#'  per protein observation.
#' @export
#'
#' @examples
#' # Single protein - Proteome discoverer format
#' sp_report_file <-
#'   system.file("extdata", "ATIC_Proteins_PD.txt", package = "tidyTPP")
#' config_file <-
#'   system.file("extdata", "ATIC_config.csv", package = "tidyTPP")
#'
#' # Import data
#' sp_tbl <- import_proteomediscoverer(sp_report_file, config_file)
#'
#' sp_tbl

import_proteomediscoverer <- function(datafile,
                                      config,
                                      path = NULL,
                                      silent = FALSE,
                                      ...){

  # Regex function to extract experiment ID:
  # digits contained in "...F(x).." e.g. "...F2.." -> 2
  experiment_id_func = \(x) sub("^.*\\.{2}F(\\d+)\\.{2}.*$", "\\1", x)


  # Give presets to custom import function
  import_custom_format(datafile, config, path,
                       table_format = "wide",
                       protein_id_col_name = "Gene.Symbol",
                       quantity_pattern = "Abundance",
                       experiment_id_func = experiment_id_func,
                       pep_n_col_name = "X..Peptides",
                       match_n_col_name = "X..PSMs",
                       silent = silent,
                       ...)

}
