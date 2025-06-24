test_that("4 protein hits by Tm are correct", {
  x <- analysis_data_4prot

  hits <-
    get_TPP_hits(
      x[c(1:3, 6:ncol(x))],
      hit_criteria = list(pvalue_max_threshold = 0.2,
                          pvalue_min_threshold = 0.15,
                          DTm_same_sign = TRUE,
                          DTm_gt_Dcontrol = TRUE,
                          slope_threshold = -0.06),
      to_export = NULL,
      to_plot = FALSE,
      silent = TRUE)

  expect_equal(
    hits,
    tibble::tibble(
      Protein_ID = "Protein_A",
      Condition = "Treated",
      Comparison = "Treated_vs_Control",
      max_adj_pvalue = 0.194889525,
      mean_melt_point = 50.9013787,
      mean_control_melt_point = 49.4093434,
      mean_diff_melt_point = 1.49203528,
      mean_control_diff_melt_point = 0.448699044
    )
  )
})

test_that("Default hits (by NPARC) are correct", {
  x <- analysis_data_4prot

  hits <-
    get_TPP_hits(
      x,
      to_export = NULL,
      to_plot = FALSE,
      silent = TRUE)

  expect_equal(
    hits,
    tibble::tibble(
      Protein_ID = "Protein_A",
      Condition = "Treated",
      Comparison = "Treated_vs_Control",
      F_scaled = 4.51772872,
      p_adj_NPARC = 0.001793434,
      max_adj_pvalue = 0.194889525,
      mean_melt_point = 50.9013787,
      mean_control_melt_point = 49.4093434,
      mean_diff_melt_point = 1.49203528,
      mean_control_diff_melt_point = 0.448699044
    ),
    tolerance = 1*10^(-6)
  )
})

test_that("No-hits behaviour is correct", {
  no_hits_data <- analysis_data_4prot[1:4, c(1:2, 6:ncol(analysis_data_4prot))]
  no_hits_data$adj_pvalue <- c(NA, NA, 1, 1)
  hits <-
    get_TPP_hits(
      no_hits_data,
      to_export = NULL,
      to_plot = FALSE,
      silent = TRUE)

  expect_equal(
    hits,
    tibble::tibble(Protein_ID = character(),
               Condition = character(),
               a = double(),
               b = double(),
               melt_point = double(),
               infl_point = double(),
               slope = double(),
               plateau = double(),
               R_sq = double(),
               Comparison = character(),
               diff_melt_point = double(),
               min_comparison_slope = double(),
               adj_pvalue = double(),
               Replicate = character(),
               mean_control_melt_point = double(),
               abs_diff_melt_control = double(),
               max_adj_pvalue = double())
  )
})
