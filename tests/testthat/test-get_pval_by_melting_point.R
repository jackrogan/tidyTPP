test_that("get_pval_by_melting_point() works as expected", {
  x <- analysis_data_4prot

  pval_x <- get_pval_by_melting_point(x[1:18])

  expect_equal(pval_x$adj_pvalue, x$adj_pvalue)
})
