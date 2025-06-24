test_that("all file types save and load identically", {

  load_Rda_and_return <- function(path){
    x <- NULL
    load(path)
    x
  }

  save_and_load <- function(x, format, func){
    file <- tempfile(fileext = format)
    export_TPP(x, file, silent = TRUE)
    y <- do.call(func, list(file))
    unlink(file)
    if(inherits(y$Replicate, "integer")){
      y$Replicate <- sprintf("%02d", y$Replicate)
    }
    tibble::as_tibble(y)
  }

  four_prot <- tibble::as_tibble(analysis_data_4prot)
  if(requireNamespace("readxl", quietly=TRUE)){
    xlsx_test <- save_and_load(four_prot, ".xlsx", readxl::read_xlsx)
    expect_equal(four_prot, xlsx_test)
  }
  csv_test <- save_and_load(four_prot, ".csv", read.csv)
  expect_equal(four_prot, csv_test)
  tsv_test <- save_and_load(four_prot, ".tsv", read.delim)
  expect_equal(four_prot, tsv_test)
  rda_test <- save_and_load(four_prot, ".Rda", load_Rda_and_return)
  expect_equal(four_prot, rda_test)
})
