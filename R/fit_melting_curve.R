# Fit melting curve to data points
# TODO - document, function to run on list (/ separate list?), error handler
fit_melting_curve <- function(data,
                              x_column = "Temp",
                              y_column = "rel_quantity",
                              protein_num = NULL,
                              protein_total = NULL,
                              silent = FALSE){
  if(any(is.null(protein_num), is.null(protein_total))) silent <- TRUE

  pre_model_data <- data.frame(x = data[[x_column]], y = data[[y_column]])
  # Fit with nls.multstart
  fit <-
    nls.multstart::nls_multstart(
      get_sigmoid_TPPTR_formula(0),
      iter = c(5,5,5),
      data = pre_model_data,
      start_lower = c(pl = -0.5, a = 0.00001, b = 0.00001),
      start_upper = c(pl = 1.5, a = 5000, b = 250),
      lower = c(pl = -0.5, a = 0.00001, b = 0.00001),
      supp_errors = "Y")

  # Get additional parameters: R2, melting point inflection point, slope
  # R2
  SS_total <-
    sum((pre_model_data$y - mean(pre_model_data$y, na.rm=TRUE))^2, na.rm=TRUE)
  SS_resid <-
    sum((pre_model_data$y - predict(fit))^2 , na.rm=TRUE)
  R_sq <- 1 - SS_resid / SS_total
  # Tm
  pars <- fit$m$getPars()
  a <- pars["a"]
  b <- pars["b"]
  pl <- pars["pl"]
  melt_point <-  a / (b - log((1 - pl) / (0.5 - pl) - 1))
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
    as.formula("y ~ ((1 - pl) / (1 + exp(b - (a / x)))) + pl")
  } else if(der == 1){
    as.formula(
      "y ~ -((1 - pl) * (exp(-(a/x - b)) * (a/x^2))/(1 + exp(-(a/x - b)))^2)"
    )
  } else if(der == 2){
    as.formula(
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
    try(uniroot(function(x, a, b, pl) eval(parse(text = formula)),
                a = a, b = b, pl = pl,
                interval = x_interval,
                tol = 0.0001),
                silent = FALSE)

  if (class(form_root) == "try-error"){
    return(NA)
  } else {
    return(form_root$root)
  }
}
