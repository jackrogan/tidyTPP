test_that("import from proteome discoverer works", {
  expect_equal(
    import_proteomediscoverer(
      system.file("extdata", "Protein_A_PD.txt", package = "tidyTPP"),
      system.file("extdata", "4_protein_config.csv", package = "tidyTPP"),
      silent = TRUE),
    tibble::tibble(quan_data_4prot[quan_data_4prot$Protein_ID == "Protein_A",]))
})
