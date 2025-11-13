x <- normalise_TPP(quan_data_20prot, silent = TRUE)

test_that("NPARC pvalue is as expected", {
  NPARC_x <-
    get_NPARC_pval(x,
                   max_cores = 1,
                   quantity_column = "rel_quantity",
                   all_stats = TRUE,
                   silent = TRUE)
  expect_equal(signif(NPARC_x[[3,10]], 2), 0.0012)
})

test_that("spline pvalue is as expected", {
  NPARC_x <-
    get_NPARC_pval(x,
                   max_cores = 1,
                   quantity_column = "rel_quantity",
                   NPARC_fit_method = "splines",
                   all_stats = TRUE,
                   silent = TRUE)
  expect_equal(signif(NPARC_x[[9,10]], 2), 0.00018)
})


test_that("NPARC with no quantity behaves as expected",{
  NPARC_x <-
    get_NPARC_pval(x[1:6],
                   max_cores = 1,
                   quantity_column = "rel_quantity",
                   silent = TRUE)
  expect_equal(NPARC_x, x[1:6])
})

test_that("NPARC distribution plots save appropriately",{
  plot_temp <- tempfile(fileext = ".png")
  expect_no_error(
    get_NPARC_pval(x,
                   max_cores = 1,
                   quantity_column = "rel_quantity",
                   NPARC_fit_method = "splines",
                   to_save = plot_temp,
                   silent = TRUE)
  )
  f_score_file <- sub("\\.png", "_Treated_Control_F_score.png", plot_temp)
  p_value_file <- sub("\\.png", "_Treated_Control_p_value.png", plot_temp)
  expect_gte(file.size(f_score_file), 30000)
  expect_gte(file.size(p_value_file), 20000)
  unlink(c(f_score_file, p_value_file))
})
