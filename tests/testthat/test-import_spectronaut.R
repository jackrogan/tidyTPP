test_that("import from spectronaut works", {
  expect_equal(
    import_spectronaut(
      system.file("extdata", "4_protein_peptide_report.csv", package = "tidyTPP"),
      system.file("extdata", "4_protein_config.csv", package = "tidyTPP"),
      silent = TRUE),
    tibble::tibble(quan_data_4prot))
})
