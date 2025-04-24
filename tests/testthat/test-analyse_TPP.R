test_that("Comparison dataframe is generated correctly for ATIC data", {
  x <- MP_data_ATIC
  reps <- unique(x$Replicate)
  conds <- unique(x$Condition)
  comparisons <- create_comparisons_tbl(conds, reps, "Control")
  expected_comparisons <- data.frame(
    Condition_01 = c("Treated", "Treated"),
    Replicate_01 = c("01", "02"),
    Condition_02 = c("Control", "Control") ,
    Replicate_02 = c("01", "02")
  )
  expect_identical(comparisons, expected_comparisons)
})

test_that("Correct melting points are calculated for ATIC treated/control", {
  x <- MP_data_ATIC
  reps <- unique(x$Replicate)
  conds <- unique(x$Condition)
  comparisons <- create_comparisons_tbl(conds, reps, "Control")
  stat_x <- find_melting_point_diffs(x, comparisons)

  expect_equal(stat_x$diff_melt_point, MP_stat_data_ATIC$diff_melt_point)
})

test_that("Correct Minimum R2 is calculated for ATIC treated/control", {
  x <- MP_stat_data_ATIC
  stat_x <- find_exp_stats(x, "min", "R_sq", "Protein_ID")

  expect_equal(stat_x$min_R_sq, x$min_R_sq)
})

test_that("Correct p-values are calculated for ATIC treated/control", {
  x <- MP_stat_data_ATIC[3:4,]

  pval_x <- calculate_binwise_pvalue(x)
  pval_x$adj_pvalue <- stats::p.adjust(pval_x$pvalue, "BH")

  expect_equal(pval_x$adj_pvalue, x$adj_pvalue)
})
