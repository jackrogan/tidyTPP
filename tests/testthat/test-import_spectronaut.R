test_that("import from spectronaut works", {
  expect_equal(
    import_spectronaut(
      system.file("extdata", "ATIC_Peptide_Report.csv", package = "tidyTPP"),
      system.file("extdata", "ATIC_config.csv", package = "tidyTPP"),
      silent = TRUE),
    tibble::tibble(quan_data_ATIC))
})
