test_that("Melting point plots with curves are generated without errors", {
  plot_return <-
    plot_melt(quan_data_4prot, analysis_data_4prot, annotate = "both")

  expect_equal(plot_return, quan_data_4prot)
})

test_that("Melting point plots without curves are generated without errors", {
  plot_return <-
    plot_melt(quan_data_4prot)

  expect_equal(plot_return, quan_data_4prot)
})

test_that("Melting plot ggplot additions work", {
  plot_return <-
    plot_melt(quan_data_4prot, analysis_data_4prot, annotate = "both",
              to_add_to_ggplot = list(ggplot2::labs(title = "ATIC Test Data")))

  expect_equal(plot_return, quan_data_4prot)
})

test_that("test ggplots are the right size", {
  temp_ggsave_file <- tempfile(fileext = ".png")
  x <- quan_data_4prot[quan_data_4prot$Protein_ID == "Protein_A",]
  y <- analysis_data_4prot[analysis_data_4prot$Protein_ID == "Protein_A",]
  plot_return <-
    suppressMessages(
       plot_melt(
        x,
        y,
        annotate = "both",
        to_add_to_ggplot =
          list(ggplot2::labs(title = "Protein_A",
                             subtitle = "Four Protein Test Data")),
        to_save = temp_ggsave_file,
        file_suffix_column = "Protein_ID"
      )
    )
  temp_ggsave_file <- sub("\\.png", "_Protein_A.png", temp_ggsave_file)
  expect_gte(file.size(temp_ggsave_file), 100000)
  unlink(temp_ggsave_file)
})
