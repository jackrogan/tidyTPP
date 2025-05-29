test_that("NPARC pvalue is as expected", {
  x <- quan_data_ATIC
  NPARC_x <-
    get_NPARC_pval(x,
                   max_cores = 1,
                   quantity_column = "rel_quantity",
                   silent = TRUE)
  expect_equal(signif(NPARC_x[[1,4]],5), 0.000031918)
})
