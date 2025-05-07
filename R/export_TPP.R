#' Export full TPP results data
#'
#' @description
#' `export_TPP()` exports a results data frame in the chosen format.
#'
#' @param TPP_data Data to export in a `data.frame` or `tibble`
#' @param file_name Character. File name to be written to. If `format` is not
#'  set, file extension will be used to decide format.
#' @param path Character. Optionally, the directory path to the location to be
#' saved in, if not the working directory.
#' @param format File save format:
#'  * \emph{xlsx:} Excel spreadsheet, using [writexl::write_xlsx()]
#'  * \emph{csv:} Comma-separated text file
#'  * \emph{tsv:} Tab-separated text file
#'  * \emph{Rda:} Rdata saved object - will loaded into R session as
#'   "TPP_results"
#'
#'  If not given, will default to "xlsx"
#'
#' @param include_columns Character vector. Column names to included in the
#'  exported data. If not given, the all columns in `TPP_data` will be exported
#' @param silent Boolean. If true, will run silently, without console output.
#'
#' @return `TPP_data`, in the format supplied
#' @export
#'
#' @examples
#' # Minimal data - analysed two-protein melt curve
#' x <- MP_stat_data_ATIC
#'
#' # Default behaviour - save a .xlsx spreadsheet of results data.
#' file <- tempfile()
#' export_TPP(x, file)
#' unlink(paste0(file, ".xlsx"))
#'
#' # Passing a filename with a .tsv extension saves as tab-separated data
#' file <- tempfile(fileext = ".tsv")
#' export_TPP(x, file)
#' unlink(file)
#'
#' # Explicitly save as .csv text file
#' file <- tempfile()
#' export_TPP(
#'   x,
#'   file,
#'   format = "csv"
#' )
#' unlink(paste0(file, ".csv"))
#'
#' # Save Rdata object with no console messages
#' file <- tempfile()
#' export_TPP(
#'   x,
#'   file,
#'   format = "Rda",
#'   silent = TRUE
#'   )
#' unlink(paste0(file, ".Rda"))
export_TPP <- function(TPP_data,
                       file_name,
                       path = NULL,
                       format = c("xlsx", "csv", "tsv", "Rda"),
                       include_columns = NULL,
                       silent = FALSE){

  TPP_results <- TPP_data

  if(!silent){
    cat("--------------------\n")
    cat("TPP Export\n")
    cat("--------------------\n")
  }

  # Choose file type by 1) file extension, 2) xlsx by default
  if(is.null(format)) format <- c("xlsx", "csv", "tsv", "Rda")
  if(identical(format, c("xlsx", "csv", "tsv", "Rda"))){
    format_pattern <- paste(format, collapse = "|")
    file_extension <-
      sub(pattern = paste0(".*\\.((", format_pattern, "))"),
          replacement = "\\1",
          x = file_name,
          ignore.case = TRUE)
    if(file_extension != file_name) format = file_extension
  }
  if(identical(format, c("xlsx", "csv", "tsv", "Rda"))){
    format <- "xlsx"
  }

  if(!silent) cat("Export format:", format, "\n")

  if(!is.null(include_columns)){
    TPP_results <- TPP_results[colnames(TPP_results) %in% include_columns]
  }

  # Combine path if given separately
  export_file <- file_name
  path <- sub("([^/])$", "\\1/", path)
  if(!is.null(path)) export_file <- paste0(path, export_file)
  if(!grepl(paste0("\\.", format, "$"), export_file)){
    export_file <- paste0(export_file, ".", format)
  }
  export_file <- sub("\\.rda$", ".Rda", export_file, ignore.case = TRUE)

  if(!silent) cat("Saving", export_file, "...\n")

  # If save method is nor recognised, don't save
  unsaved <- c(
    "Export format not recognised, please use 'xlsx', 'csv', 'tsv' or 'Rda'.\n",
    "Not saved.\n")
  save_func <- \(x, y) cat(unsaved, sep = "")
  if(toupper(format) == "RDA") save_func <- \(x, y) save(x, file = y)
  if(toupper(format) == "CSV") {
    save_func <- \(x, y) utils::write.table(x, y, sep = ",")
  }
  if(toupper(format) == "TSV"){
    save_func <- \(x, y) utils::write.table(x, y, sep = "\t")
  }

  if(toupper(format) == "XLSX") save_func <- writexl::write_xlsx

  do.call(save_func, args = list(TPP_results, export_file))

  if(!silent) cat("Saved.\n--------------------\n")

  TPP_results
}
