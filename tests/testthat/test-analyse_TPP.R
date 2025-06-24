test_that("create_comparisons_tbl() works", {
  x <- analysis_data_4prot[analysis_data_4prot$Protein_ID == "Protein_A",]
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

test_that("create_control_comparison_tbl() works", {
  x <- analysis_data_4prot[analysis_data_4prot$Protein_ID == "Protein_A",]
  reps <- unique(x$Replicate)
  conds <- unique(x$Condition)
  controls <- create_control_comparison_tbl(reps, "Control")
  expected_controls <- data.frame(
    Condition_01 = c("Control"),
    Replicate_01 = c("01"),
    Condition_02 = c("Control") ,
    Replicate_02 = c( "02")
  )
  expect_identical(controls, expected_controls)
})

test_that("find_melting_point_diffs() works", {
  x <- analysis_data_4prot[analysis_data_4prot$Protein_ID == "Protein_A",]
  reps <- unique(x$Replicate)
  conds <- unique(x$Condition)
  comparisons <- create_comparisons_tbl(conds, reps, "Control")
  controls <- create_control_comparison_tbl(reps, "Control")
  comparisons <- rbind(comparisons, controls)
  stat_x <- find_melting_point_diffs(x, comparisons)

  expect_equal(stat_x$diff_melt_point, x$diff_melt_point)
})

test_that("find_exp_stats() works for R2", {
  x <- analysis_data_4prot[analysis_data_4prot$Protein_ID == "Protein_A",]
  stat_x <- find_exp_stats(x, "min", "R_sq", "Protein_ID")

  expect_equal(stat_x$min_R_sq, x$min_R_sq)
})
