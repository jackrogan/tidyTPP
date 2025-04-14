#' Fit Melting Curve for Experimental Groups
#'
#' @description
#' `fit_melt_by_experiment()` fits sigmoidal melting curves to all the proteins
#' in `data` using [fit_melting_curve()], splitting into individual curves
#' identified with supplied columns.
#'
#' Both serial and parallel processing (using `PSOCK`) of the individual curves
#' are supported.
#'
#' @inheritParams fit_melting_curve
#' @param experiment_cols A character vector of column names, which will be used
#' to split the data into constituent experimental groups. This should include
#' a column that uniquely identifies proteins in the data
#' @param max_cores Integer. The maximum number of cores to parallelise the
#' curve fitting over, if less than the maximum number of cores available.
#' If `max_cores` is 1, then curve fitting will be run serially in
#' a single process
#' @param ... Further arguments to be passed to [fit_melting_curve()] and
#' [nls_multstart()]
#'
#' @return A `tibble` with 1 row and \emph{7 + n} variables (where
#' \emph{n} is the length of `experiment_cols`):
#' * ID columns from `experiment_cols`
#' * `a`: Parameter `a` from fitted sigmoidal curve.
#' * `b`: Parameter `b` from fitted sigmoidal curve.
#' * `plateau`: Plateau from fitted sigmoidal curve.
#' * `melt_point`:  Melting point (\eqn{T_m}), temperature at which 50 % of the
#'    protein is destabilised
#' * `infl_point`: Inflection point (\eqn{T_{infl}}), temperature at which the
#'    second derivative \eqn{f''(T) = 0}
#' * `slope`: Sigmoidal curve slope, \emph{i.e.} \eqn{f'(T_{infl})}
#' * `R_sq`: \eqn{R^2} for the fitted sigmoidal curve
#' @export
#'
#' @examples
#' # Minimal data - ATIC melt curve
#' x <- quan_data_ATIC
#'
#' # Fit melting curves - 2-core processing
#' fit_melt_by_experiment(x, max_cores = 2)
#'
#' # Single-core (serial) processing
#' fit_melt_by_experiment(x, max_cores = 1)
#'
#' # Specify experiment details to separate.
#' fit_melt_by_experiment(
#'   x,
#'   experiment_cols = c("Protein_ID", "Condition", "Replicate"),
#'   max_cores = 2
#' )
#'
fit_melt_by_experiment <-
  function(data,
           experiment_cols = c("Protein_ID", "Condition", "Replicate"),
           silent = FALSE,
           max_cores = 8,
           ...){

    # Get dots to pass to cluster
    dots <- list(...)

    # Supplying null values to split should result in not splitting
    if(length(experiment_cols) != 0){
      expr_cols <- paste("~", paste(experiment_cols, collapse = " + "))
      data_sep <- split(data, stats::as.formula(expr_cols))
    } else {
      data_sep <- list(data)
    }

    # Fitting function
    fit_curve_x <- function(x){
      params_x <-
        do.call(fit_melting_curve,
                c(list("data" = data_sep[[x]],
                       "protein_num" = x,
                       "protein_total" = total,
                       "silent" = silent),
                       dots))
      if(!is.null(params_x)) {
        params_x <- cbind(data_sep[[x]][experiment_cols][1,], params_x)
      }
      params_x
    }

    # Run fitting function (TODO in parallel?)
    total = length(data_sep)
    if(!silent) cat("\nFitting", total, "proteins...\n")
    time_total <- proc.time()

    # Run in parallel
    if(max_cores > 1){
      n_cores <- min(parallel::detectCores(), max_cores)

      # Per operation assume 0.2 s
      if(!silent) {
        cat("Running on", n_cores, "cores.\nEstimated total process time:",
            round(0.65 + total * 0.23 / n_cores, 2), "s\n")
      }
      cl <- parallel::makeCluster(n_cores)
      parallel::clusterExport(cl, list(
        "fit_melting_curve",
        "get_sigmoid_TPPTR_formula",
        "get_sigmoid_formula_root",
        "experiment_cols",
        "silent",
        "dots",
        "data_sep",
        "total"
      ),
      envir = environment())

      time_total <- proc.time()
      params <-
        parallel::parLapply(cl, 1:length(data_sep), fit_curve_x) |>
        Reduce(rbind, x = _)

      parallel::stopCluster(cl)
    # Run serially
    } else {
      if(!silent) {
        cat("Estimated total process time:", round(total*0.23, 2), "s\n")
      }
      params <-
        lapply(1:length(data_sep), fit_curve_x) |>
        Reduce(rbind, x = _)
    }
    if(!silent) cat(nrow(params), "of", total, "fitted successfully.\n")

    time_out <- proc.time() - time_total
    cat("\nTotal elapsed time:", time_out[["elapsed"]], "s\n")

    tibble::as_tibble(params)
  }
