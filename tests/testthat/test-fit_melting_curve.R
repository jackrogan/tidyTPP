test_that("ATIC curve fit returns tibble", {
  x_sep <- split(quan_data_ATIC, ~ Condition + Replicate)
  total = length(x_sep)
  params <-
    lapply(1:total,
           function(i) {
             fit_melting_curve(x_sep[[i]], protein_num = i, protein_total = total)
           }) |>
    Reduce(rbind, x = _)

  expect_s3_class(params, "tbl_df")
})
