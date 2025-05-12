test_that("ATIC hits by Tm are correct", {
  hits <-
    get_TPP_hits(
      MP_stat_data_ATIC,
      hit_criteria = list(DTm_same_sign = TRUE,
                          DTm_gt_Dcontrol = TRUE),
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
