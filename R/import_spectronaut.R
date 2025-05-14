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
import_spectronaut <- function(datafile,
                               config,
                               path = NULL,
                               silent = FALSE,
                               ...){

  # Regex function to extract experiment ID:
  # digits following "Sample_"
  experiment_id_func = \(x) gsub(".*Sample_(\\d+)_.*", "\\1", x)

  # Regext function to extract canonical sequence from precursor ID:
  # remove leading and trailing "_", charge (as ".x") where x is a digit,
  # and modification labels in "[]"
  seq_col_func = \(x) gsub("_|(\\.\\d+)|(\\[.*\\])", "", x)

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
