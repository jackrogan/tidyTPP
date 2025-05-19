# Function to get statistics experiment-wise
find_exp_stats <- function(x_tbl, stat_func, stat_column, experiment_cols){

  exp_formula <-
    stats::as.formula(
      paste("cbind(", paste(stat_column, collapse = ", "), ") ~" ,
            paste(experiment_cols, collapse = " + ")))

  for(f in stat_func){
    f_tbl <- aggregate(exp_formula, x_tbl, FUN = f)
    stat_col_names <- colnames(f_tbl) %in% stat_column
    colnames(f_tbl)[stat_col_names] <- paste(f, stat_column, sep = "_")
    x_tbl <- merge(x_tbl, f_tbl)
  }
  x_tbl
}

# Additional functions
min_abs <- function(x) min(abs(x))
all_same_sign <- function(x) all(x > 0) | all(x < 0)
