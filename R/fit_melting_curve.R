#' Fit Sigmoidal Melting Curve and Find Parameters
#'
#' @description
#' `fit_melting_curve()` fits observed protein quantity data to a sigmoidal
#' melting point curve, given by:
#' \deqn{  f(T)=\frac{1-plateau}{1+e^{-(\frac{a}{T}-b)}}+plateau }
#'
#' Curve fitting uses [nls_multstart] to estimate a series of starting
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
#' @param ... Further arguments to be passed to [nls_multstart()]
#'
#' @return A `tibble` with 1 row and 7 variables:
#' * `a`: Parameter `a` from fitted sigmoidal curve.
#' * `b`: Parameter `b` from fitted sigmoidal curve.
#' * `plateau`: Plateau from fitted sigmoidal curve.
#' * `melt_point`:  Melting point (\eqn{T_m}), temperature at which 50 % of the
#'    protein is destabilised
#' * `infl_point`: Inflection point (\eqn{T_{infl}}), temperature at which the
#'    second derivative \eqn{f''(T) = 0}
#' * `slope`: Sigmoidal curve slope, \emph{i.e.} \eqn{f'(T_{infl})}
#' * `R_sq`: \eqn{R^2} for the fitted sigmoidal curve
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
#' # Minimal data - ATIC melt curve
#' x <- quan_data_ATIC
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
                              ...){
  if(any(is.null(protein_num), is.null(protein_total))) silent <- TRUE

  pre_model_data <- data.frame(x = data[[x_column]], y = data[[y_column]])

  # Allow dots to override default arguments for nls
  nls_dots <- list(...)
  nls_defaults <- list("formula" = get_sigmoid_TPPTR_formula(0),
                       "iter" = c(5,5,5),
                       "data" = pre_model_data,
                       "start_lower" = c(pl = -0.5, a = 0.00001, b = 0.00001),
                       "start_upper" = c(pl = 1.5, a = 5000, b = 250),
                       "lower" = c(pl = -0.5, a = 0.00001, b = 0.00001),
                       "supp_errors" = "Y")
  nls_args <- list()
  for(i in 1:length(nls_defaults)){
    if(!names(nls_defaults)[i] %in% names(nls_dots)) {
      nls_args <- c(nls_args, nls_defaults[i])
    }
  }
  nls_args <- c(nls_args, nls_dots)

  # Fit with nls.multstart
  fit <-
    try(
      do.call(nls.multstart::nls_multstart, nls_args),
      silent = show_errors
    )
  if(inherits(fit, "try-error")){
    if(!silent){
      if(protein_num == protein_total) cat("\n")
    }
    return(NULL)
  }

  # Get additional parameters: R2, melting point inflection point, slope
  # R2
  SS_total <-
    sum((pre_model_data$y - mean(pre_model_data$y, na.rm=TRUE))^2, na.rm=TRUE)
  SS_resid <-
    sum((pre_model_data$y - stats::predict(fit))^2 , na.rm=TRUE)
  R_sq <- 1 - SS_resid / SS_total
  # Tm
  pars <- fit$m$getPars()
  a <- pars["a"]
  b <- pars["b"]
  pl <- pars["pl"]
  melt_point <-  suppressWarnings(a / (b - log((1 - pl) / (0.5 - pl) - 1)))
  if(is.na(melt_point)) melt_point <- NA
  # Tinfl
  T_interval <- c(min(pre_model_data$x), max(pre_model_data$x))
  infl_point <-
    get_sigmoid_formula_root(get_sigmoid_TPPTR_formula(2),
                             a = a, b = b, pl = pl,
                             x_interval = T_interval)
  # Slope
  x <- infl_point
  if(!is.na(infl_point)){
    slope <- eval(parse(text = get_sigmoid_TPPTR_formula(1)))
  } else {
    slope <- NA
  }

  if(!silent) {
    progress = floor(10 * protein_num / protein_total)
    cat("\r|", strrep("=", progress), strrep(" ", 10 - progress), "| ",
        protein_num, " of ", protein_total, sep = "")
    if(protein_num == protein_total) cat("\n")
  }

  tibble::tibble(a, b, plateau = pl, melt_point, infl_point, slope, R_sq)
}

# Function to return melting curve sigmoid formula and derivatives
get_sigmoid_TPPTR_formula <- function(der = 0){
  if(der == 0){
    stats::as.formula("y ~ ((1 - pl) / (1 + exp(b - (a / x)))) + pl")
  } else if(der == 1){
    stats::as.formula(
      "y ~ -((1 - pl) * (exp(-(a/x - b)) * (a/x^2))/(1 + exp(-(a/x - b)))^2)"
    )
  } else if(der == 2){
    stats::as.formula(
      "y ~ -((1 - pl) * 1 * (exp(-(a/x - b)) * (a/x^2) *
      (a/x^2) - exp(-(a/x - b)) * (a * (2 * x)/(x^2)^2)) /
      (1 + exp(-(a/x - b)))^2 - (1 - pl) * 1 * (exp(-(a/x - b)) *
      (a/x^2)) * (2 * (exp(-(a/x - b)) * (a/x^2) *
      (1 + exp(-(a/x - b)))))/((1 + exp(-(a/x - b)))^2)^2)"
    )
  }
}

# Function to get root of supplied formula
get_sigmoid_formula_root <- function(formula, a, b, pl, x_interval){
  form_root <-
    try(stats::uniroot(function(x, a, b, pl) eval(parse(text = formula)),
                       a = a, b = b, pl = pl,
                       interval = x_interval,
                       tol = 0.0001),
        silent = TRUE)

  if(inherits(form_root, "try-error")){
    return(NA)
  } else {
    return(form_root$root)
  }
}
