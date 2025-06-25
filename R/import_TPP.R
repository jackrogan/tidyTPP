#' Import TPP data from chosen proteomics tabular data format into tidier tibble
#' format
#'
#' @description
#'  `import_TPP()` imports tabular proteomics data from a delimited
#'  text file or \emph{.xlsx} spreadsheet and prepares it for TPP data analysis,
#'  by choosing between functions that transform specific formats.
#'
#'  * Columns specified using the given arguments are selected from the data,
#'    according to the vendor file format
#'
#'  * Experiment condition combinations (treatment, temperature, \emph{etc.}) are
#'    loaded from a separate config file and matched to recorded data
#'
#'  * The necessary transformations are performed so the data can be analysed for
#'    \emph{thermal protein profiling (TPP)} data: combining data points so that
#'    a single data point is presented for each experiment combination and
#'    recording the number of peptide hits for each protein as far as that's
#'    possible from the data - this allows the use of minimum matches as a
#'    quality control.
#'
#'  * The data is reshaped into a tidier long format `tibble`, with one
#'    observation per row, one column per statistic.
#'
#' @inheritParams import_custom_format
#' @param import_format Character. Specify a format to import with an
#'  \emph{import_} function, or custom specification:
#'
#'  * "Spectronaut" or "SN" - peptide report from \emph{Spectronaut} DIA protein
#'  quantification (\emph{.csv}): [import_spectronaut]
#'
#'  * "ProteomeDiscover" or "PD" - exported proteins table from
#'  \emph{Thermo Proteome Discoverer} DDA quantification (\emph{.txt}):
#'  [import_proteomediscoverer]
#'
#'  * "Custom" - custom specification: [import_custom_format]
#'
#' @param ... Arguments to be passed to `import_` functions, and then to
#' `read.table()` or `read_excel()`
#'
#' @return
#'  A tibble giving the TPP protein quantity data relative to the lowest
#'  temperature, \eqn{T_1}, as \emph{rel_quantity}, in a tidier format: 1 row
#'  per protein observation.
#' @export
#'
#' @examples
#' # Four-protein example - Spectronaut format
#' sp_report_file <-
#'   system.file("extdata", "4_protein_peptide_report.csv", package = "tidyTPP")
#' config_file <-
#'   system.file("extdata", "4_protein_config.csv", package = "tidyTPP")
#'
#' # Import data
#' sp_tbl <-
#'   import_TPP(sp_report_file, config_file, import_format = "Spectronaut")
#'
#' sp_tbl
import_TPP <- function(datafile,
                       config,
                       path = NULL,
                       import_format = "Spectronaut",
                       silent = FALSE,
                       ...){

  null_import_func <- function(datafile, config, silent, ...){
    if(!silent){
      cat("\n--------------------\n")
      cat("TPP Data Import\n")
      cat("--------------------\n")
      cat("Invalid import format. Please supply one of the following:\n")
      cat(" - Spectronaut (or SN)",
          " - ProteomeDiscoverer (or PD)",
          " - Custom", sep = "\n")
      cat("--------------------\n")
    }
    return(NULL)
  }

  import_func <-
    switch(toupper(import_format),
           "SPECTRONAUT" = import_spectronaut,
           "SN" = import_spectronaut,
           "PROTEOMEDISCOVERER" = import_proteomediscoverer,
           "PD" = import_proteomediscoverer,
           "CUSTOM" = import_custom_format,
           null_import_func
    )

  do.call(import_func,
          list(datafile = datafile,
               config = config,
               path = path,
               silent = silent,
               ...))
}
