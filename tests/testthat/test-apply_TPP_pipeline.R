test_that("pipeline runs without errors", {

  temp_loc <- tempdir()
  source_report <-
    system.file("extdata", "normP_report.csv", package = "tidyTPP")
  source_config <-
    system.file("extdata", "ATIC_config.csv", package = "tidyTPP")

  file.copy(c(source_report, source_config), temp_loc)

  expect_no_error(
    apply_TPP_pipeline(
      paste0("normP_report.csv"),
      paste0("ATIC_config.csv"),
      path = temp_loc,
      to_plot = TRUE,
      max_cores = 1
    )
  )

  unlink(
    c(
      paste0(temp_loc, "\\", "normP_report.csv"),
      paste0(temp_loc, "\\", "ATIC_config.csv"),
      paste0(temp_loc, "\\", "normP_report_Results.xlsx"),
      paste0(temp_loc, "\\", "normP_report_Hits.xlsx"))
  )
})
