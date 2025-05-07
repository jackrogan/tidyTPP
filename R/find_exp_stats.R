# Function to get statistics experiment-wise
find_exp_stats <- function(x_tbl, stat_func, stat_column, experiment_cols){
  apply_stat <- function(x){
    for(f in stat_func){
      for(i in stat_column){
        x[paste(f, i, sep = "_")] <- do.call(get(f), list(x[, i]))
      }
    }
    x
  }

  exp_formula <-
    stats::as.formula(paste("~" , paste(experiment_cols, collapse = " + ")))
  x_tbl <-
    split(x_tbl, exp_formula) |>
    lapply(apply_stat) |>
    Reduce(rbind, x = _)

  x_tbl
}

# Additional functions
min_abs <- function(x) min(abs(x))
all_same_sign <- function(x) all(x > 0) | all(x < 0)
