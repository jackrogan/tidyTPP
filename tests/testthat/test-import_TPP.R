test_that("spectronaut is chosen correctly", {
  expect_equal(
    import_TPP(
      system.file("extdata", "4_protein_peptide_report.csv", package = "tidyTPP"),
      system.file("extdata", "4_protein_config.csv", package = "tidyTPP"),
      import_format = "spectronaut",
      silent = TRUE),
    tibble::tibble(quan_data_4prot))
})

test_that("Proteome Dsicoverer is chosen correctly", {
  expect_equal(
    import_TPP(
      system.file("extdata", "Protein_A_PD.txt", package = "tidyTPP"),
      system.file("extdata", "4_protein_config.csv", package = "tidyTPP"),
      import_format = "PrOtEomEDiScovERer",
      silent = TRUE),
    tibble::tibble(quan_data_4prot[quan_data_4prot$Protein_ID == "Protein_A",]))
})

test_that("Custom import is chosen correctly", {
  PD_prot_A_data <-
    system.file("extdata", "Protein_A_PD.txt", package = "tidyTPP")
  config <-
    system.file("extdata", "4_protein_config.csv", package = "tidyTPP")

  experiment_id_func = \(x) sub("^.*\\.{2}F(\\d+)\\.{2}.*$", "\\1", x)

  expect_equal(
    import_TPP(PD_prot_A_data,
               config,
               import_format = "CUSTOM",
               table_format = "wide",
               protein_id_col_name = "Gene.Symbol",
               quantity_pattern = "Abundance",
               experiment_id_func = experiment_id_func,
               pep_n_col_name = "X..Peptides",
               match_n_col_name = "X..PSMs",
               silent = TRUE),
    tibble::tibble(quan_data_4prot[quan_data_4prot$Protein_ID == "Protein_A",])
  )
})

test_that("Unrecognized format is handled correctly", {
  expect_null(import_TPP("x", "y", import_format = "no_format", silent = FALSE))
})
