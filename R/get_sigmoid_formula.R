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
