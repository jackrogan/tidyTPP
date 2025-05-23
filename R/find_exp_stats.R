# Function to get statistics experiment-wise
find_exp_stats <- function(x_tbl, stat_func, stat_column, experiment_cols){

  exp_formula <-
    stats::as.formula(
      paste("cbind(", paste(stat_column, collapse = ", "), ") ~" ,
            paste(experiment_cols, collapse = " + ")))

  for(i in 1:length(stat_func)){
    f <- stat_func[[i]]
    f_name <- names(stat_func)[i]
    if(is.null(f_name)){
      if(inherits(f, "character")) {
        f_name <- f
      } else {
        f_name <- paste0("f", i)
      }
    }
    f_tbl <- stats::aggregate(exp_formula, x_tbl, FUN = f)
    stat_col_names <- colnames(f_tbl) %in% stat_column

    colnames(f_tbl)[stat_col_names] <-
      paste(f_name, stat_column, sep = "_")
    x_tbl <- merge(x_tbl, f_tbl)
  }
  x_tbl
}
