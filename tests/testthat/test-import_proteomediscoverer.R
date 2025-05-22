test_that("import from spectronaut works", {
  expect_equal(
    import_proteomediscoverer(
      system.file("extdata", "ATIC_Proteins_PD.txt", package = "tidyTPP"),
      system.file("extdata", "ATIC_config.csv", package = "tidyTPP"),
      silent = TRUE),
    tibble::tibble(quan_data_ATIC))
})
