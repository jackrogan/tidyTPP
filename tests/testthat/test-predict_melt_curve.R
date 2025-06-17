test_that("curve is predicted from nls result", {
  x <- quan_data_ATIC
  x <- x[x$Condition == "Control" & x$Replicate == "02",]
  x <- mask_column(x, "Temp", "x")
  x <- mask_column(x, "rel_quantity", "y")
  x_formula <- get_sigmoid_TPPTR_formula(0)
  x_predic <-
    stats::nls(x_formula,
               x,
               start =  c(pl = 0.005, a = 1500, b = 10),
               algorithm = "port")
  x_curve <- predict_melt_curve(x_predic, 10, 37, 67)

  expect_equal(nrow(x_curve), 10)
})

test_that("curve is predicted from list of nls result", {
  x <- quan_data_ATIC
  x <- mask_column(x, "Temp", "x")
  x <- mask_column(x, "rel_quantity", "y")
  x_list <- split(x, ~ Condition + Replicate)
  x_formula <- get_sigmoid_TPPTR_formula(0)

  get_x_predicted <- function(x) {
    x_predic <-
      stats::nls(x_formula,
                x,
                start =  c(pl = 0.005, a = 1500, b = 10),
                algorithm = "port")
  }
  x_predic_list <- lapply(x_list, get_x_predicted)

  x_curve <- predict_melt_curve(x_predic_list, 10, 37, 67)

  expect_equal(nrow(x_curve), 40)
})
