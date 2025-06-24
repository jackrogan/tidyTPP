test_that("column masking functions as expected", {
  x <- quan_data_4prot
  x$quantity = 1
  y <- mask_column(x, "rel_quantity", "quantity")
  z <- mask_column(y, "quantity", "rel_quantity")
  expect_equal(x, z)
  y_names <- colnames(x)
  y_names[y_names == "quantity"] <- "quantity_renamed"
  y_names[y_names == "rel_quantity"] <- "quantity"
  expect_equal(y_names, colnames(y))
})
