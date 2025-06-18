test_that("NPARC pvalue is as expected", {
  x <- quan_data_ATIC
  NPARC_x <-
    get_NPARC_pval(x,
                   max_cores = 1,
                   quantity_column = "rel_quantity",
                   silent = TRUE)
  expect_equal(signif(NPARC_x[[1,4]],5), 0.000031918)
})

test_that("NPARC with no quantity behaves as expected",{
  x <- quan_data_ATIC[1:6]
  NPARC_x <-
    get_NPARC_pval(x,
                   max_cores = 1,
                   quantity_column = "rel_quantity",
                   silent = TRUE)
  expect_equal(as.data.frame(NPARC_x), x)
})

test_that("NPARC distribution plots save appropriately",{
  plot_temp <- tempfile(fileext = ".png")
  x <- quan_data_ATIC
  expect_no_error(
    get_NPARC_pval(x,
                   max_cores = 1,
                   quantity_column = "rel_quantity",
                   to_save = plot_temp,
                   silent = TRUE)
  )
  f_score_file <- sub("\\.png", "_Treated_Control_F_score.png", plot_temp)
  p_value_file <- sub("\\.png", "_Treated_Control_p_value.png", plot_temp)
  expect_equal(file.size(f_score_file) >= 60000, TRUE)
  expect_equal(file.size(p_value_file) >= 50000, TRUE)
  unlink(c(f_score_file, p_value_file))
})
