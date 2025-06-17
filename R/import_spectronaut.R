#' Import TPP data from Spectronaut report into tidier tibble format.
#'
#' @inheritParams import_custom_format
#' @param datafile Character. Specify a Spectronaut .csv report filename
#'  containing the protein quantity data.
#'
#' @return A tibble giving the TPP protein quantity data relative to the lowest
#'  temperature, \eqn{T_1}, as \emph{rel_quantity}, in a tidier format: 1 row
#'  per protein observation.
#' @export
#'
#' @examples
#'
#' # Single protein - Spectronaut format
#' sp_report_file <-
#'   system.file("extdata", "ATIC_Peptide_Report.csv", package = "tidyTPP")
#' config_file <-
#'   system.file("extdata", "ATIC_config.csv", package = "tidyTPP")
#'
#' # Import data
#' sp_tbl <- import_spectronaut(sp_report_file, config_file)
#'
#' sp_tbl

import_spectronaut <- function(datafile,
                               config,
                               path = NULL,
                               silent = FALSE,
                               ...){

  # Regex function to extract experiment ID:
  # digits following "Sample_"
  experiment_id_func <- \(x) gsub(".*Sample_(\\d+)_.*", "\\1", x)

  # Regex function to extract canonical sequence from precursor ID:
  # remove leading and trailing "_", charge (as ".x") where x is a digit,
  # and modification labels in "[]"
  seq_col_func <- \(x) gsub("_|(\\.\\d+)|(\\[.*\\])", "", x)

  # Give presets to custom import function
  import_custom_format(datafile, config, path,
                       table_format = "wide",
                       protein_id_col_name = "PG.Genes",
                       quantity_pattern = "PG.Quantity",
                       experiment_id_func = experiment_id_func,
                       seq_col_name = "EG.PrecursorId",
                       seq_col_func = seq_col_func,
                       silent = silent,
                       ...)

}
