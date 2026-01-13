#' Fit Sigmoidal Melting Curve and Find Parameters
#'
#' @description
#' `fit_melting_curve()` fits observed protein quantity data to a sigmoidal
#' melting point curve, given by:
#' \deqn{  f(T)=\frac{1-plateau}{1+e^{-(\frac{a}{T}-b)}}+plateau }
#'
#' Curve fitting iterates calls to [nls] until the starting parameters result in
#' convergence to a least-squares estimate for best fit.
#' Arguments given to `nls()` by default are:
#' * `start = c(pl = 0, a = 550, b = 10)`
#' * `lower = c(pl = 0, a = 0.00001, b = 0.00001)`
#' * `upper = c(pl = 1.5, a = 15000, b = 250)`
#' * `na.action = na.omit`
#' * `algorithm = "port"`
#' * `control = stats::nls.control(maxiter=50))`
#'
#' Alternatively, [nls_multstart] can be used to estimate a series of starting
#' parameters for the non-linear model using a gridstart approach.
#' Arguments given to `nls_multstart()` by default are:
#' * `iter = c(5,5,5)`
#' * `start_lower = c(pl = -0.5, a = 0.00001, b = 0.00001)`
#' * `start_upper = c(pl = 1.5, a = 5000, b = 250)`
#' * `lower = c(pl = -0.5, a = 0.00001, b = 0.00001)`
#'
#' \eqn{R^2} is calculated using a pseudo-\eqn{R^2} approach for the non-linear
#' model, \eqn{T_m} from solving \eqn{f(T) = 0.5}, \eqn{T_{infl}} from solving
#' \eqn{f''(T) = 0}, and slope by evaluating \eqn{f'(T_{infl})}.
#'
#' @param data A data frame (or [tibble]) containing proteomics data from a
#'  thermal protein profiling (TPP) experiment. This must contain the columns:
#'  * Experiment temperature in \eqn{\degree C}
#'  * Protein quantity measurement corresponding to the non-denatured fraction
#'
#' @param x_column Character. Name for the column containing x-axis
#' (temperature) values in `data`
#' @param y_column Character. Name for the column containing y-axis
#' (protein quantity) values in `data`
#' @param protein_num,protein_total Numeric values to be displayed as part of
#' the progress bar if necessary, ignored if silent = `TRUE`
#'
#'  \emph{e.g.} for `num` = 30, `total` = 256:
#'
#'  `> |=         | 30 of 256`
#'
#' @param silent Boolean. If `TRUE`, no console output is shown
#' @param show_errors Boolean. If `TRUE`, errors in calling [nls_multstart()]
#'  will be shown in console output
#' @param curve_fit_method Character. Method for sigmoid curve fitting -
#'  either \emph{"nls"} for a repeated [nls()] approach (default) or
#'  \emph{"nls.multstart"} for [nls_multstart()].
#' @param ... Further arguments to be passed to [nls()] or [nls_multstart()]
#'
#' @return A `tibble` with 1 row and 13 variables:
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
#' * `RSS`: Residual sum of squares (\emph{RSS})
#' * `sigma`: Residual standard deviation (\eqn{\sigma})
#' * `n_coeffs`: Number of coefficients fitted
#' * `n_obs`: Number of observations in curve
#' * `log_lik`: log-likelihood (\eqn{\ell})
#' * `AICc`: Corrected Akaike information criterion (\emph{AICc})
#'
#' @references
#'  Schellman J. A., The thermodynamics of solvent exchange
#'  \emph{Biopolymers} 34, 1015-1026 (1994)
#'
#'  Savitski M. M. \emph{et al.}, Tracking cancer drugs in living cells by
#'  thermal profiling of the proteome. \emph{Science}, 346: 1255784 (2014)
#' @export
#'
#' @examples
#' # Minimal data - 4-protein melt curve (Protein A only)
#' x <- quan_data_4prot[quan_data_4prot$Protein_ID == "Protein_A",]
#'
#' # Separate condition and replicate
#' x_sep <- split(x, ~ Condition + Replicate)
#'
#' # Fit melting curves
#' total = length(x_sep)
#' params <-
#'   lapply(1:total,
#'     function(i) {
#'       fit_melting_curve(x_sep[[i]], protein_num = i, protein_total = total)
#'     }) |>
#'   Reduce(rbind, x = _)
#'
#' params
#'
fit_melting_curve <- function(data,
                              x_column = "Temp",
                              y_column = "rel_quantity",
                              protein_num = NULL,
                              protein_total = NULL,
                              silent = FALSE,
                              show_errors = FALSE,
                              curve_fit_method = c("nls", "nls.multstart"),
                              ...){

  if(any(is.null(protein_num), is.null(protein_total))) silent <- TRUE

  pre_model_data <- data.frame(x = data[[x_column]], y = data[[y_column]])

  if("nls" %in% curve_fit_method){

    # Allow dots to override default arguments for nls
    nls_dots <- list(...)
    nls_defaults <- list("formula" = get_sigmoid_TPPTR_formula(0),
                         "data" = pre_model_data,
                         "start" = c(pl = 0, a = 550, b = 10),
                         "max_attempts" = 100,
                         "lower" = c(pl = 0, a = 0.00001, b = 0.00001),
                         "upper" = c(pl = 1.5, a = 15000, b = 250),
                         "na.action" = stats::na.omit,
                         "algorithm" = "port",
                         "control" = stats::nls.control(maxiter=50))
    nls_args <-
      nls_defaults[sapply(names(nls_defaults), \(x) !x %in% names(nls_dots))]
    nls_args <- c(nls_args, nls_dots)
    nls_args <- nls_args[names(nls_args) %in% names(formals(stats::nls))]
    # Use try and repeat

    fit <- do.call(repeat_nls, nls_args)
  } else if("nls.multstart" %in% curve_fit_method){
    # Allow dots to override default arguments for nls.multstart
    nls_dots <- list(...)
    nls_defaults <- list("formula" = get_sigmoid_TPPTR_formula(0),
                         "iter" = c(5,5,5),
                         "data" = pre_model_data,
                         "start_lower" = c(pl = -0.5, a = 0.00001, b = 0.00001),
                         "start_upper" = c(pl = 1.5, a = 15000, b = 250),
                         "lower" = c(pl = -0.5, a = 0.00001, b = 0.00001),
                         "supp_errors" = "Y")
    nls_args <-
      nls_defaults[sapply(names(nls_defaults), \(x) !x %in% names(nls_dots))]
    nls_args <- c(nls_args, nls_dots)
    nls_args <-
      nls_args[
        names(nls_args) %in% names(formals(nls.multstart::nls_multstart))]

    # Fit with nls.multstart
    fit <-
      try(
        do.call(nls.multstart::nls_multstart, nls_args),
        silent = !show_errors
      )

  } else {
    if(!silent) cat("\nNo curve fit method selected - cannot fit.")
    return(NULL)
  }

  if(inherits(fit, "try-error")){
    if(!silent){
      if(protein_num == protein_total) cat("\n")
    }
    return(NULL)
  }


  if(!silent) {
    progress = floor(10 * protein_num / protein_total)
    cat("\r|", strrep("=", progress), strrep(" ", 10 - progress), "| ",
        protein_num, " of ", protein_total, sep = "")
    if(protein_num == protein_total) cat("\n")
  }

  tibble::tibble(get_model_fit_stats(fit))
}

# Function to attempt nls fit and repeat
repeat_nls <- function(start = c(pl = 0, a = 550, b = 10),
                       max_attempts = 100,
                       show_errors = FALSE,
                       ...){

  i <- 0
  repeat_fit <- TRUE

  while (repeat_fit){
    mod_start <- start * (1 + stats::runif(1, -0.5, 0.5))
    fit <-
      try(do.call(stats::nls, list("start" = mod_start, ...)),
          silent = !show_errors)
    i <- i + 1
    repeat_fit <- inherits(fit, "try-error") & i < max_attempts
  }

  fit
}
