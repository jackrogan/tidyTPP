#' Parameters for ATIC Protein Melting Curve
#'
#' A dataset containing the parameters generated from fitting a melting curve
#' to raw  proteomics data (DIA) for the protein ATIC:
#'
#' @format A data frame with 4 rows and 10 variables:
#' \describe{
#'   \item{Protein_ID}{Gene name identifying the protein}
#'   \item{Condition}{Treatment given to this sample set}
#'   \item{Replicate}{Replicate number for this sample set}
#'   \item{`a`}{Parameter `a` from sigmoidal curve fitting}
#'   \item{`b`}{Parameter `b` from sigmoidal curve fitting}
#'   \item{melt_point}{Melting point, calculated as the model temperature
#'    resulting in 50 % of the protein destabilised}
#'   \item{infl_point}{Inflection point, calculated as the model temperature
#'    for which \eqn{f''(T_{infl}) = 0 }}
#'   \item{slope}{Sigmoidal curve slope, calculated as \eqn{f'(T_{infl}) }}
#'   \item{plateau}{The plateau of the modeled sigmoidal curve}
#'   \item{R_sq}{\eqn{R^2} for the fitted sigmoidal curve}
#' }
MP_data_ATIC <-
  data.frame(
    "Protein_ID" = rep("ATIC", times = 4),
    "Condition" = c("Control", "Control", "Treated", "Treated"),
    "Replicate" = c("01", "02", "01", "02"),
    "a" = c(967.148, 1101.623, 1267.360, 1180.769),
    "b" =  c(19.72895, 22.24584, 24.60575, 23.08822),
    "melt_point" = c(49.02177, 49.52040, 51.51124, 51.14162),
    "infl_point" = c(48.52645, 49.12545, 51.17007, 50.76260),
    "slope" = c(-0.1016437, -0.1132118, -0.1200861, -0.1137089),
    "plateau" = c(0.000000000, 0.000000000, 0.001091275, 0.000000000),
    "R_sq" = c(0.9972189, 0.9941574, 0.9993804, 0.9981549)
  )

#' Quantity data for ATIC Protein Melting Curve
#'
#' A dataset containing the parameters generated from fitting a melting curve
#' to raw  proteomics data (DIA) for the protein ATIC:
#'
#' @format A data frame with 10 rows and 5 variables:
#' \describe{
#'   \item{Protein_ID}{Gene name identifying the protein}
#'   \item{Pep_N}{Number of unique peptides identified for protein}
#'   \item{Match_N}{Number of unique peptide spectral matches (PSMs)
#'    found for protein}
#'   \item{Condition}{Treatment given to this sample set}
#'   \item{Replicate}{Replicate number for this sample set}
#'   \item{Temp}{Experiment temperature in \eqn{\degree C}}
#'   \item{rel_quantity}{Quantity of protein measured, relative to \eqn{T_1}}
#'   \item{raw_quantity}{Raw quantity of protein measured}
#' }
quan_data_ATIC <-
  data.frame(
    "Protein_ID" = rep("ATIC", times = 40),
    "Pep_N" = rep(as.integer(40), times = 40),
    "Match_N" = rep(as.integer(62), times = 40),
    "Condition" = rep(c("Control", "Treated"), each = 20),
    "Replicate" = rep(c("01","02"), each = 10, times = 2),
    "Temp" = rep(c(37, 41, 44, 47, 50, 53, 56, 59, 63, 67), times = 4),
    "rel_quantity" = c(1.000000000, 0.914750055, 0.895068033, 0.689273176,
                       0.415987281, 0.194960365, 0.045819291, 0.016755450,
                       0.025820331, 0.009859520, 1.000000000, 1.021613752,
                       0.922667249, 0.700220756, 0.439630669, 0.195675589,
                       0.052509818, 0.014621917, 0.010663716, 0.007950098,
                       1.000000000, 1.257181252, 1.329592830, 1.276262079,
                       0.996100563, 0.492405634, 0.139852090, 0.031413770,
                       0.010516553, 0.017891194, 1.000000000, 0.966318210,
                       0.990162739, 0.842307100, 0.661046993, 0.393178545,
                       0.115087087, 0.027304649, 0.010410181, 0.008438892),
    "raw_quantity" = c(7913465.50, 7238843.00, 7083090.00, 5454539.50,
                       3291901.00, 1542812.12,  362589.38,  132593.67,
                        204328.30,   78022.97, 7737944.00, 7905190.00,
                       7139547.50, 5418269.00, 3401837.50, 1514126.75,
                        406318.03,  113143.58,   82515.23,   61517.41,
                       5856743.00, 7362987.50, 7787083.50, 7474739.00,
                       5833905.00, 2883893.25,  819077.75,  183982.38,
                         61592.75,  104784.12, 7475048.00, 7223275.00,
                       7401514.00, 6296286.00, 4941358.00, 2939028.50,
                        860281.50,  204103.56,   77816.60,   63081.12)
  )

# Minimal list of nls models: ATIC Melting curve - control 01
model_ATIC_1 <-
  quan_data_ATIC[quan_data_ATIC$Condition == "Control" &
                   quan_data_ATIC$Replicate == "01",] |>
  nls(rel_quantity ~ ((1 - plateau) / (1 + exp(b - (a / Temp)))) + plateau,
      data = _,
      start = list(plateau = 0, a = 600, b = 10))
