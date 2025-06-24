test_that("incompatible custom specs are handled: 1. No Quantities", {

  PD_prot_A_data <-
    system.file("extdata", "Protein_A_PD.txt", package = "tidyTPP")
  config <-
    system.file("extdata", "4_protein_config.csv", package = "tidyTPP")

  experiment_id_func = \(x) sub("^.*\\.{2}F(\\d+)\\.{2}.*$", "\\1", x)

  expect_null(
    import_custom_format(PD_prot_A_data, config,
                         table_format = "wide",
                         protein_id_col_name = "Gene.Symbol",
                         quantity_pattern = "THIS_NAME_LEFT_WRONG",
                         experiment_id_func = experiment_id_func,
                         pep_n_col_name = "X..Peptides",
                         match_n_col_name = "X..PSMs",
                         silent = TRUE)
  )

})

test_that("incompatible custom specs are handled: 2. No Protein IDs", {

  PD_prot_A_data <-
    system.file("extdata", "Protein_A_PD.txt", package = "tidyTPP")
  config <-
    system.file("extdata", "4_protein_config.csv", package = "tidyTPP")

  experiment_id_func = \(x) sub("^.*\\.{2}F(\\d+)\\.{2}.*$", "\\1", x)

  expect_null(
    import_custom_format(PD_prot_A_data, config,
                         table_format = "wide",
                         protein_id_col_name = "THIS_NAME_LEFT_WRONG",
                         quantity_pattern = "Abundance",
                         experiment_id_func = experiment_id_func,
                         pep_n_col_name = "X..Peptides",
                         match_n_col_name = "X..PSMs",
                         silent = TRUE)
  )

})

test_that("incompatible custom specs are handled: 3. No Experiment IDs", {

  PD_prot_A_data <-
    system.file("extdata", "Protein_A_PD.txt", package = "tidyTPP")
  config <-
    system.file("extdata", "4_protein_config.csv", package = "tidyTPP")

  experiment_id_func = \(x) sub("^.*\\.{2}F(\\d+)\\.{2}.*$", "\\1", x)

  expect_null(
    import_custom_format(PD_prot_A_data, config,
                         table_format = "long",
                         protein_id_col_name = "Gene.Symbol",
                         quantity_pattern = "Abundance",
                         experiment_col_name = "THIS_NAME_LEFT_WRONG",
                         pep_n_col_name = "X..Peptides",
                         match_n_col_name = "X..PSMs",
                         silent = TRUE)
  )

})
