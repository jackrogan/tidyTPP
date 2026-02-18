#' Import TPP data from proteomics tabular data into tidier tibble format.
#'
#' @description
#' `import_custom_format()` imports tabular proteomics data from a delimited
#' text file or \emph{.xlsx} spreadsheet and prepares it for TPP data analysis,
#' allowing a specification for a given proteomics software file type to be
#' easily and consistently produced:
#'
#' * Columns specified using the given arguments are selected from the data,
#'   according to the vendor file format
#'
#' * Experiment condition combinations (treatment, temperature, \emph{etc.}) are
#'   loaded from a separate config file and matched to recorded data
#'
#' * The necessary transformations are performed so the data can be analysed for
#'   \emph{thermal protein profiling (TPP)} data: combining data points so that
#'   a single data point is presented for each experiment combination and
#'   recording the number of peptide hits for each protein as far as that's
#'   possible from the data - this allows the use of minimum matches as a
#'   quality control.
#'
#' * The data is reshaped into a tidier long format `tibble`, with one
#'   observation per row, one column per statistic.
#'
#' @param datafile Character. Specify a filename for a delimited text file or
#' \emph{.xlsx} spreadsheet containing the protein quantity data.
#' @param config  Either a character, specifying the file name of a delimited
#'  file with experiment configuration, to be opened with
#'  [read.table] / [readxl::read_excel] or a `data.frame` / `tibble` with
#'  configuration information.
#'
#'  This should take the form of a table with the column names
#'  \emph{Experiment}, \emph{Condition}, \emph{Replicate} and
#'  \emph{Temp}.
#' @param path Character. Optionally, the path to a directory containing the
#'  data and configuration files. If given, this will be appended to both paths.
#' @param table_format One of "wide" or "long". Are the quantity data in the
#'  file to import presented in a wide format - one column per temperature point
#'  observation, or long - one quantity column, one experimental identifier
#'  column, \emph{e.g.}
#'
#'  \emph{Wide:}
#'
#'  | Quantity_Sample_01.raw | Quantity_Sample_02.raw | Quantity_Sample_03.raw |
#'  | ----------------------:| ----------------------:| ----------------------:|
#'  | 5434942                | 5394587                | 4523683                |
#'
#'  \emph{Long:}
#'
#'  | Experiment    | Quantity |
#'  | -------------:| --------:|
#'  | Sample_01.raw | 5434942  |
#'  | Sample_02.raw | 5394587  |
#'  | Sample_03.raw | 4523683  |
#'
#' @param quantity_pattern Character. Regular expression that identifies protein
#'  quantity column(s)
#' @param protein_id_col_name Character. The exact name of a column in the data
#' from which a unique protein ID will be taken.
#' @param experiment_col_name Character. If `table_format = "long"`, the exact
#'  name of a column in the data from which experiment ID should be taken.
#' @param quantity_col_name Character. If `table_format = "wide"`, the name of
#'  a column in the data from which protein quantity values should be taken.
#'  Will override `quantity_pattern` if both are supplied
#' @param experiment_id_func Optional function to transform vector of column
#'  names for quantity columns into experiment IDs matching the contents of the
#'  `Experiment` column from the `config` data. If not given, samples are
#'  assigned a number based on order of appearance in the data.
#' @param pep_n_col_name Optional character. The exact name of a column in the
#'  data from which number of peptide matches should be taken. If given, assumes
#'  single record per protein and does not use sequences.
#' @param match_n_col_name Optional character. The exact name of a column in the
#'  data from which number of PSM matches should be taken. Used if present and
#'  `pep_col_name` is given.
#' @param seq_col_name Optional character. The exact name of a column in the
#'  data from which peptide sequences should be taken to calculate number of
#'  peptide matches. Will be ignored if `pep_col_name` is given. If not given,
#'  number of peptide IDs per protein will be estimated from number of identical
#'  rows in data per protein.
#' @param seq_col_func Optional function to transform peptide sequence column to
#'  extract peptide canonical sequence (used to differentiate peptide matches
#'  per protein from total spectrum matches)
#' @param silent Boolean. If true, will run silently, without console output.
#' @param ... Arguments to pass to `read.table()` or `read_excel()`
#'
#' @return A tibble giving the TPP protein quantity data relative to the lowest
#'  temperature, \eqn{T_1}, as \emph{rel_quantity}, in a tidier format: 1 row
#'  per protein observation.
#' @export
#'
#' @examples
#' # Four proteins - Spectronaut format
#' sp_report_file <-
#'   system.file("extdata", "4_protein_peptide_Report.csv", package = "tidyTPP")
#' config_file <-
#'   system.file("extdata", "4_protein_config.csv", package = "tidyTPP")
#'
#' # Regex function to extract experiment ID:
#' # digits following "Sample_"
#' experiment_id_func = \(x) gsub(".*Sample_(\\d+)_.*", "\\1", x)
#'
#' # Regex function to extract canonical sequence from precursor ID:
#' # remove leading and trailing "_", charge (as ".x") where x is a digit,
#' # and modification labels in "[]"
#' seq_col_func = \(x) gsub("_|(\\.\\d+)|(\\[.*\\])", "", x)
#'
#' # Load data using paramaters for exported spectronaut report .csv
#' sp_tbl <- import_custom_format(
#'   # Report file
#'   datafile = sp_report_file,
#'   # Configuration
#'   config = config_file,
#'   # Spectronaut uses wider data format
#'   table_format = "wide",
#'   # Use quantity column pattern for multiple quantity columns
#'   protein_id_col_name = "PG.Genes",
#'   quantity_pattern = "PG.Quantity",
#'   experiment_id_func = experiment_id_func,
#'   # Use sequence column to generate number of peptide matches
#'   seq_col_name = "EG.PrecursorId",
#'   seq_col_func = seq_col_func
#' )
#'
#' sp_tbl

import_custom_format <- function(datafile,
                                 config,
                                 path = NULL,
                                 table_format = c("wide", "long"),
                                 protein_id_col_name = NULL,
                                 quantity_pattern = NULL,
                                 experiment_col_name = NULL,
                                 quantity_col_name = NULL,
                                 experiment_id_func = NULL,
                                 pep_n_col_name = NULL,
                                 match_n_col_name = NULL,
                                 seq_col_name = NULL,
                                 seq_col_func = NULL,
                                 silent = FALSE,
                                 ...){

  if(!silent){
    cat("\n--------------------\n")
    cat("TPP Data Import\n")
    cat("--------------------\n")
    cat("Read in:\n")
    if(is.data.frame(datafile)){
      cat("Data table supplied\n")
    } else{
      cat(datafile, "\n")
    }
    if(is.data.frame(config)){
      cat("Configuration table supplied\n")
    } else{
      cat(config, "\n")
    }
  }

  # Data/Config import
  data_in_tbl <- multiformat_import(datafile, path)
  config_tbl <- multiformat_import(config, path)

  if(!silent) cat("--------------------\n")

  # Function to get colnames that match arguments only if not null
  colname_grep <- function(pattern, return_type = "v"){
    if(!is.null(pattern)) {
      if(return_type == "l") {
        grepl(pattern, colnames(data_in_tbl))
      } else {
        grep(pattern, colnames(data_in_tbl), value = TRUE)
      }
    } else {
      NULL
    }
  }

  # Default wide format
  table_format <- table_format[1]

  quan_cols <- NULL
  # Required columns:
  # Protein ID, Peptide Precursor (if present), PG raw quantities
  # if(!is.null(quantity_pattern)){
  #   quan_cols <- grep(quantity_pattern, colnames(data_in_tbl), value = TRUE)
  # }
  quan_cols <- colname_grep(quantity_pattern)
  # if(table_format == "long" & !is.null(quantity_col_name)){
  #   quan_cols <- quantity_col_name
  # }
  if(table_format == "long") quan_cols <- colname_grep(quantity_col_name)
  if(is.null(quan_cols) | length(quan_cols) == 0){
    if(!silent) cat("No quantity columns identified!\n")
    return(NULL)
  }
  # protein_id_col_name <-
  #   colnames(data_in_tbl)[grepl(protein_id_col_name, colnames(data_in_tbl))]
  protein_id_col_name <- colname_grep(protein_id_col_name)
  if(length(protein_id_col_name) == 0){
    if(!silent) cat("No protein ID columns identified!\n")
    return(NULL)
  }
  # if(!is.null(experiment_col_name)){
  #   experiment_col_name <-
  #     colnames(data_in_tbl)[grepl(experiment_col_name, colnames(data_in_tbl))]
  # }
  experiment_col_name <- colname_grep(experiment_col_name)
  if(table_format == "long" & length(experiment_col_name) == 0){
    if(!silent) cat("No experiment ID columns identified!\n")
    return(NULL)
  }
  if(all(!is.null(seq_col_name),
         !any(colname_grep(seq_col_name, return_type = "l")),
         !silent)){
    cat("Sequence column given does not match a column in data!\n")
    seq_col_name <- NULL
  }

  pep_n_col_name <- colname_grep(pep_n_col_name)
  match_n_col_name <- colname_grep(match_n_col_name)
  TPP_cols <- c(protein_id_col_name,
                pep_n_col_name,
                match_n_col_name,
                seq_col_name,
                experiment_col_name,
                quan_cols)

  TPP_tbl <- data_in_tbl[TPP_cols]
  colnames(TPP_tbl)[1] <- "Protein_ID"

  # If Peptide matches and PSM matches are given, use those, otherwise assume
  # multiple rows per protein and aggregate
  if(!is.null(pep_n_col_name)){
    colnames(TPP_tbl)[colnames(TPP_tbl) == pep_n_col_name] <- "Pep_N"
    if(!is.null(match_n_col_name)){
      colnames(TPP_tbl)[colnames(TPP_tbl) == match_n_col_name] <- "Match_N"
    }

  } else {

    if(!is.null(seq_col_name)) colnames(TPP_tbl)[2] <- "Sequence"

    if(!is.null(seq_col_func) & !is.null(seq_col_name)) {
      TPP_tbl$Sequence <- seq_col_func(TPP_tbl$Sequence)
    }


    # If peptide sequences are given, get number of peptides and number of total
    # matches. Else use total matches per protein

    pep_N_tbl <- TPP_tbl
    if(!is.null(seq_col_name)){
      pep_N_tbl$Match_N <- 0
      match_N_tbl <-
        stats::aggregate(Match_N ~ Protein_ID + Sequence ,
                         pep_N_tbl, FUN = length)
      match_N_tbl$Pep_N <- 1
      pep_N_tbl <-
        stats::aggregate(cbind(Pep_N, Match_N) ~ Protein_ID,
                         match_N_tbl, FUN = sum)
    } else{
      pep_N_tbl$Pep_N <- 0
      pep_N_tbl <- stats::aggregate(Pep_N ~ Protein_ID, pep_N_tbl, FUN = length)
    }

    # Reform data table
    TPP_tbl <- TPP_tbl[c("Protein_ID", experiment_col_name, quan_cols)]
    TPP_tbl <- merge(pep_N_tbl, TPP_tbl)
  }
  # Reduce to single protein value
  TPP_tbl <- stats::aggregate(. ~ Protein_ID, TPP_tbl, FUN = \(x) x[1])
  # Drop missing protein values
  TPP_tbl <- TPP_tbl[TPP_tbl$Protein_ID != "",]


  # Transform to long format if not already
  if(table_format == "wide"){
    if(!silent) cat("Pivoting to long table...\n")
    TPP_tbl <-
      stats::reshape(TPP_tbl,
                     direction = "long",
                     varying = quan_cols,
                     timevar = "Experiment",
                     v.names = "raw_quantity",
                     idvar = "Protein_ID",
                     times = quan_cols)

  }

  # Transform Experiment column
  if(!is.null(experiment_id_func)){
    if(!silent) cat("Transforming experiment names...\n")
    # TPP_tbl$Experiment <- as.numeric(experiment_id_func(TPP_tbl$Experiment))
    TPP_tbl$Experiment <- experiment_id_func(TPP_tbl$Experiment)
  }

  # Match to Experiment values
  if(!silent) cat("Matching to experiment config data...\n")
  TPP_tbl <- merge(config_tbl, TPP_tbl)
  # Get T1 values. relative quantity
  if(!silent) cat("Finding relative quantity values...\n")
  # TODO allow user-selected reference channel label
  if("Pooled" %in% TPP_tbl$Temp){
    if(!silent) cat("Using 'Pooled' reference channel label.\n")
    TPP_tbl <- get_quantity_relative_to_reference(TPP_tbl)
  } else {
    TPP_tbl$rel_quantity <- TPP_tbl$raw_quantity
  }

  TPP_tbl <- get_quantity_relative_to_T1(TPP_tbl)
  TPP_tbl$Temp <- as.numeric(TPP_tbl$Temp)

  TPP_tbl <-
    TPP_tbl[order(TPP_tbl$Protein_ID,
                  TPP_tbl$Condition,
                  TPP_tbl$Replicate,
                  TPP_tbl$Temp),]
  # Replicate should be discrete value
  TPP_tbl$Replicate <- sprintf("%02d", TPP_tbl$Replicate)

  cols_out <- c("Protein_ID", "Pep_N", "Match_N",
                "Condition", "Replicate", "Temp",
                "rel_quantity", "raw_quantity")
  cols_out <- cols_out[cols_out %in% colnames(TPP_tbl)]
  TPP_tbl <- TPP_tbl[cols_out]
  n_proteins <- length(unique(TPP_tbl$Protein_ID))

  if(!silent) {
    cat("Data imported.\n")
    cat("Found", n_proteins, "proteins.\n")
    cat("--------------------\n")
  }

  tibble::as_tibble(TPP_tbl)
}

get_quantity_relative_to_T1 <- function(TPP_tbl){
  T1_tbl <- TPP_tbl[order(TPP_tbl$Temp),]
  T1_tbl <- stats::aggregate(rel_quantity ~ Protein_ID + Condition + Replicate,
                             T1_tbl, FUN = \(x) x[1])
  colnames(T1_tbl)[4] <- "T1_quantity"
  TPP_tbl <- merge(TPP_tbl, T1_tbl)
  TPP_tbl$rel_quantity <- TPP_tbl$rel_quantity / TPP_tbl$T1_quantity
  TPP_tbl
}

get_quantity_relative_to_reference <- function(TPP_tbl, ref_name = "Pooled"){
  pooled_tbl <- TPP_tbl[TPP_tbl$Temp == ref_name,]
  pooled_tbl$ref_quantity <- pooled_tbl$raw_quantity
  pooled_tbl <-
    pooled_tbl[c("Protein_ID", "Condition", "Replicate", "ref_quantity")]
  TPP_tbl <- TPP_tbl[TPP_tbl$Temp != ref_name,]
  TPP_tbl <- merge(TPP_tbl, pooled_tbl)
  TPP_tbl$rel_quantity <- TPP_tbl$raw_quantity / TPP_tbl$ref_quantity
  #TPP_tbl[colnames(TPP_tbl) != "ref_quantity"]
  TPP_tbl
}

multiformat_import <- function(file, path = NULL, ...){
  if(!is.data.frame(file)){

    # (include "/" if not given)
    if(!is.null(path)){
      path <- sub("([^/])$", "\\1/", path)
      file <- paste0(path, file)
    }
    if(file.exists(file)){

      # Choose read in method from extension
      file_extension <- toupper(sub(".*(\\.[^\\.]+$)", "\\1", file))
      read_func <-
        switch(file_extension,
               ".XLSX" = readxl::read_xlsx,
               ".XLS" = readxl::read_xls,
               ".CSV" = utils::read.csv,
               utils::read.delim
        )
      file <- do.call(read_func, c(file, list(...)))
    } else {
      file <-  NULL
    }
  }
  file
}
