test_that("Melting point plots with curves are generated without errors", {
  plot_return <-
    plot_melt(quan_data_ATIC, MP_data_ATIC, annotate = "both")

  expect_equal(plot_return, quan_data_ATIC)
})

test_that("Melting point plots without curves are generated without errors", {
  plot_return <-
    plot_melt(quan_data_ATIC)

  expect_equal(plot_return, quan_data_ATIC)
})

test_that("Melting plot ggplot additions work", {
  plot_return <-
    plot_melt(quan_data_ATIC, MP_data_ATIC, annotate = "both",
              to_add_to_ggplot = list(ggplot2::labs(title = "ATIC Test Data")))

  expect_equal(plot_return, quan_data_ATIC)
})

test_that("test ggplots are the right size", {
  temp_ggsave_file <- tempfile(fileext = ".png")
  plot_return <-
    suppressMessages(
       plot_melt(
        quan_data_ATIC,
        MP_data_ATIC,
        annotate = "both",
        to_add_to_ggplot = list(ggplot2::labs(title = "ATIC Test Data")),
        to_save = temp_ggsave_file,
        file_suffix_column = "Protein_ID"
      )
    )
  temp_ggsave_file <- sub("\\.png", "_ATIC.png", temp_ggsave_file)
  expect_true(file.size(temp_ggsave_file) >= 210000 &
                file.size(temp_ggsave_file) < 240000)
  unlink(temp_ggsave_file)
})
