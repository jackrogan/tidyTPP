test_that("ATIC hits by Tm are correct", {
  hits <-
    get_TPP_hits(
      MP_stat_data_ATIC[c(1:3, 6:ncol(MP_stat_data_ATIC))],
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
      Protein_ID = "ATIC",
      Condition = "Treated",
      Comparison = "Treated_vs_Control",
      max_adj_pvalue = 0.142924515,
      mean_melt_point = 51.3264300,
      mean_control_melt_point = 49.2710850,
      mean_diff_melt_point = 2.05534500,
      mean_control_diff_melt_point = 0.4986300
    )
  )
})

test_that("Default hits (by NPARC) are correct", {
  hits <-
    get_TPP_hits(
      MP_stat_data_ATIC,
      to_export = NULL,
      to_plot = FALSE,
      silent = TRUE)

  expect_equal(
    hits,
    tibble::tibble(
      Protein_ID = "ATIC",
      Condition = "Treated",
      Comparison = "Treated_vs_Control",
      F_scaled = 8.78931,
      p_adj_NPARC = 3.191767e-05,
      max_adj_pvalue = 0.142924515,
      mean_melt_point = 51.3264300,
      mean_control_melt_point = 49.2710850,
      mean_diff_melt_point = 2.05534500,
      mean_control_diff_melt_point = 0.4986300
    )
  )
})

test_that("No-hits behaviour is correct", {
  no_hits_data <- MP_stat_data_ATIC[c(1:2, 6:ncol(MP_stat_data_ATIC))]
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
               adj_pvalue = double(),
               Replicate = character(),
               mean_control_melt_point = double(),
               abs_diff_melt_control = double(),
               max_adj_pvalue = double())
  )
})
