test_that("pipeline runs without errors", {

  temp_loc <- tempdir()
  source_report <-
    system.file("extdata", "4_protein_peptide_report.csv", package = "tidyTPP")
  source_config <-
    system.file("extdata", "4_protein_config.csv", package = "tidyTPP")

  file.copy(c(source_report, source_config), temp_loc)

  expect_no_error(
    apply_TPP_pipeline(
      paste0("4_protein_peptide_report.csv"),
      paste0("4_protein_config.csv"),
      path = temp_loc,
      to_plot = TRUE,
      max_cores = 1
    )
  )

  unlink(
    c(
      paste0(temp_loc, "\\", "4_protein_peptide_report.csv"),
      paste0(temp_loc, "\\", "4_protein_config.csv"),
      paste0(temp_loc, "\\", "4_protein_peptide_report_Results.xlsx"),
      paste0(temp_loc, "\\", "4_protein_peptide_report_Hits.xlsx"))
  )
})
