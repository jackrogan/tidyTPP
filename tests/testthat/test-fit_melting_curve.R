test_that("4 protein curve fit returns tibble", {
  x_sep <- split(quan_data_4prot, ~ Protein_ID + Condition + Replicate)
  total = length(x_sep)
  params <-
    lapply(1:total,
           function(i) {
             fit_melting_curve(x_sep[[i]],
                               protein_num = i,
                               protein_total = total,
                               silent = TRUE)
           }) |>
    Reduce(rbind, x = _)

  expect_s3_class(params, "tbl_df")
})

test_that("4 protein curve fit returns tibble with multstart", {
  x_sep <- split(quan_data_4prot, ~ Protein_ID + Condition + Replicate)
  total = length(x_sep)
  params <-
    lapply(1:total,
           function(i) {
             fit_melting_curve(x_sep[[i]],
                               protein_num = i,
                               protein_total = total,
                               silent = TRUE,
                               curve_fit_method = "nls.multstart")
           }) |>
    Reduce(rbind, x = _)

  expect_s3_class(params, "tbl_df")
})
