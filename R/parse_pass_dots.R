# Function to deal with arguments to pass to functions:
# Pass only arguments named by function
parse_pass_dots <- function(func, add_args = list(), ...){
  dots <- list(...)
  func_dots <- dots[names(dots) %in% names(formals(func))]
  do.call(func, c(add_args, func_dots))
}
