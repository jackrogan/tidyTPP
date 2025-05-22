test_that("calculate_binwise_pvalue() works as expected", {
  x <- MP_stat_data_ATIC[3:4,]

  pval_x <- calculate_binwise_pvalue(x)
  pval_x$adj_pvalue <- stats::p.adjust(pval_x$pvalue, "BH")

  expect_equal(pval_x$adj_pvalue, x$adj_pvalue)
})
