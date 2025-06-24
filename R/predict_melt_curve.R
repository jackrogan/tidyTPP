#' Predict Protein Melting Curve
#'
#' @param fit_parameters Parameters for plotting melting curve. Must include at
#'  minimum `plateau`, `a`, `b`. Accepted formats:
#'  * A `data.frame`, column titles must include `plateau`, `a`, `b`.
#'  * A vector of numeric parameters, of length 3, with the order `plateau`, `a`, `b`.
#'  * A model, created using [nls()], from which parameters can be extracted.
#'  * A list of [nls()] models, from which parameters can be extracted and made
#'  into a `data.frame`
#' @param n_predict An integer, the number of predictions to make per melting point
#'  curve
#' @param low_T A number, the lower temperature limit for predicted curves
#' @param high_T A number, the upper temperature limit for predicted curves
#' @param T_seq Optionally, a numeric vector of temperature values to
#'  predict; overrides `n_predict`, `low_T` and `high_T`
#'
#' @return A `tibble` (from [tibble()]) including original `data.frame` columns,
#'  if any, `Temp`, and paired `model_quantity`
#' @export
#'
#' @examples
#' # 4-protein example (Protein A only)
#' x <- analysis_data_4prot[analysis_data_4prot$Protein_ID == "Protein_A",]
#' predict_melt_curve(x)
#'
#' predict_melt_curve(x, n_predict = 100, low_T = 36, high_T = 66)
predict_melt_curve <-
  function(fit_parameters,
           n_predict = 100,
           low_T = 37,
           high_T = 67,
           T_seq = NULL){

    # Determine fit parameter type, get data.frame version
    if(is.data.frame(fit_parameters)) {
      parameter_tbl <- fit_parameters

    # Vector of parameters
    } else if(is.numeric(fit_parameters) & length(fit_parameters) == 3){
      parameter_tbl <- .get_params_df_from_vector(fit_parameters)

    # Model from nls
    } else if(inherits(fit_parameters, "nls")){
      parameter_tbl <- .get_params_df_from_nls(fit_parameters)

    # List of models from nls
    } else if(is.list(fit_parameters) & inherits(fit_parameters[[1]], "nls")){
      parameter_tbl <-
        lapply(fit_parameters, .get_params_df_from_nls) |>
        Reduce(rbind, x = _)

      # Return early if none of the above
    } else {
      return()
    }

    # Generate modeled values
    n_models <- nrow(parameter_tbl)
    if(!is.null(T_seq)){
      parameter_tbl <- parameter_tbl[rep(1:n_models, each = length(T_seq)),]
      parameter_tbl$Temp <- rep(T_seq, times = n_models)
    } else{
      parameter_tbl <- parameter_tbl[rep(1:n_models, each = n_predict),]
      parameter_tbl$Temp <-
        rep(seq(low_T, high_T, length.out = n_predict), times = n_models)
    }
    parameter_tbl$model_quantity <-
      ((1 - parameter_tbl$plateau) /
         (1 + exp(parameter_tbl$b - (parameter_tbl$a / parameter_tbl$Temp)))) +
      parameter_tbl$plateau


    tibble::tibble(parameter_tbl)
  }

# Helper function to return model parameters vector as data.frame
.get_params_df_from_vector <- function(params){
  data.frame(plateau = params[1],
             a = params[2],
             b = params[3])
}

# Helper function to return parameters from nls model as data.frame
.get_params_df_from_nls <- function(model){
  model$m$getPars()[c("pl", "a", "b")] |>
    .get_params_df_from_vector()
}

