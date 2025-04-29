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

  ATIC <- tibble::as_tibble(MP_stat_data_ATIC)
  if(requireNamespace("readxl", quietly=TRUE)){
    xlsx_test <- save_and_load(ATIC, ".xlsx", readxl::read_xlsx)
    expect_equal(ATIC, xlsx_test)
  }
  csv_test <- save_and_load(ATIC, ".csv", read.csv)
  expect_equal(ATIC, csv_test)
  tsv_test <- save_and_load(ATIC, ".tsv", read.delim)
  expect_equal(ATIC, tsv_test)
  rda_test <- save_and_load(ATIC, ".Rda", load_Rda_and_return)
  expect_equal(ATIC, rda_test)
})
