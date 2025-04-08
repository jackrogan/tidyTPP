#' Import TPP data from Spectronaut report into tidier tibble format.
#'
#' @param datafile Character. Specify a Spectronaut .csv report filename
#'  containing the protein quantity data.
#' @param config  Either a character, specifying the file name of a text .csv
#'  file with experiment configuration, or literal data that can be passed to
#'  [readr::read_csv]. This should take the form of a table with the
#'  column names \emph{Experiment}, \emph{Condition}, \emph{Replicate} and
#'  \emph{Temp}.
#' @param path Character. Optionally, the path to a directory containing the
#'  data and configuration files. If given, this will be appended to both paths.
#' @param silent Boolean. If true, will run silently, without console output.
#'
#' @return A tibble giving the TPP protein quantity data relative to the lowest
#'  temperature, \eqn{T_1}, as \emph{rel_quantity}, in a tidier format: 1 row
#'  per protein observation.
#' @export
#'
#' @importFrom rlang .data
import_spectronaut <- function(datafile,
                               config,
                               path = NULL,
                               silent = FALSE){

  if(!silent){
    cat("--------------------\n")
    cat("Spectronaut Import\n")
    cat("--------------------\n")
    cat("Read in:\n")
    if(is.data.frame(datafile)){
      cat("Spectronaut report table supplied\n")
    } else{
      cat(datafile, "\n")
    }
    if(is.data.frame(config)){
      cat("Configuration table supplied\n")
    } else{
      cat(config, "\n")
    }
  }

  # If path is given, add to input names
  # (include "/" if not given)
  if(!is.null(path)){
    path <- sub("([^/])$", "\\1/", path)
    if(!is.data.frame(datafile)) datafile <- paste0(path, datafile)
    if(!is.data.frame(config)) config <- paste0(path, config)
  }

  # Check for missing files, import data
  if(!is.data.frame(datafile)){
    if(!file.exists(datafile)) {
      cat("Can't find datafile:\n")
      cat(datafile, "\n")
      return(NULL)
    }
  }
  if(!is.data.frame(config)){
    if(!file.exists(config)) {
      cat("Can't find config:\n")
      cat(config, "\n")
      return(NULL)
    }
  }
  SN_data_tbl <- readr::read_csv(datafile, show_col_types = FALSE)
  config_tbl <- readr::read_csv(config, col_types = "cccn", show_col_types = FALSE)
  # Required columns:
  # Protein Group, Gene, Peptide Precursor, PG raw quantities
  SN_data_tbl <-
    dplyr::select(SN_data_tbl,
                  'Protein_ID' = 'PG.Genes',
                  'EG.PrecursorId',
                  dplyr::contains("raw.PG.Quantity")) |>
    dplyr::distinct() |>
    # Separate Precursor sequence, charge
    tidyr::separate_wider_delim(cols = 'EG.PrecursorId',
                                delim = ".",
                                names = c("Pep.Sequence", "Pep.Charge")) |>
    # Create peptide match and PSM match count
    dplyr::group_by(.data$Protein_ID) |>
    dplyr::summarise('Pep_N' = length(unique(.data$Pep.Sequence)),
                     'Match_N' = dplyr::n(),
                     dplyr::across(dplyr::contains("raw.PG.Quantity"),
                                   ~ utils::head(.x, 1))) |>
    # Remove any data missing in any experiment
    tidyr::drop_na() |>
    # Combine config and data tibbles to give temp-quantity pairs
    tidyr::pivot_longer(cols = dplyr::contains("raw.PG.Quantity"),
                        names_to = "Experiment",
                        values_to = "raw_quantity") |>
    dplyr::mutate('Experiment' = sub("^\\[(\\d+)\\].*$", "\\1", .data$Experiment)) |>
    dplyr::full_join(config_tbl, by = "Experiment") |>
    dplyr::group_by(.data$Protein_ID, .data$Condition, .data$Replicate) |>
    # Get T1 temperature and generate relative quantity column
    dplyr::mutate('T1_quantity' = utils::head(.data$raw_quantity, 1),
                  'rel_quantity' = .data$raw_quantity / .data$T1_quantity) |>
    dplyr::ungroup()


  if(!silent){ cat("--------------------\n") }

  dplyr::select(SN_data_tbl,
                'Protein_ID', 'Pep_N', 'Match_N', 'Condition',
                'Replicate', 'Temp', 'rel_quantity', 'raw_quantity')
}
