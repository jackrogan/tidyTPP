#' Quantity data for 4 Protein Melting Curves
#'
#' A dataset containing the raw proteomics data (DIA) for 4 proteins
#'
#' @format A data frame with 160 rows and 8 variables:
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
"quan_data_4prot"

#' Statistics for 4 Protein Melting Curve
#'
#' A dataset containing the statistics and parameters generated from fitting a
#' melting curve to normalised proteomics data (DIA) for the 4 proteins and
#' comparing melting curves across measurements:
#'
#' @format A data frame with 16 rows and 19 variables:
#' \describe{
#'   \item{Protein_ID}{Gene name identifying the protein}
#'   \item{Condition}{Treatment given to this sample set}
#'   \item{Replicate}{Replicate number for this sample set}
#'   \item{F_scaled}{Scaled F-statistic from NPARC analysis}
#'   \item{p_adj_NPARC}{p-value from NPARC analysis, FDR adjusted using
#'    \emph{Benjamini-Hochberg} correction}
#'   \item{`a`}{Parameter `a` from sigmoidal curve fitting}
#'   \item{`b`}{Parameter `b` from sigmoidal curve fitting}
#'   \item{melt_point}{Melting point, calculated as the model temperature
#'    resulting in 50 % of the protein destabilised}
#'   \item{infl_point}{Inflection point, calculated as the model temperature
#'    for which \eqn{f''(T_{infl}) = 0 }}
#'   \item{slope}{Sigmoidal curve slope, calculated as \eqn{f'(T_{infl}) }}
#'   \item{plateau}{The plateau of the modeled sigmoidal curve}
#'   \item{R_sq}{\eqn{R^2} for the fitted sigmoidal curve}
#'   \item{Comparison}{Measurements compared for \eqn{\Delta T_m}}
#'   \item{diff_melt_point}{Melting point difference (\eqn{\Delta T_m})}
#'   \item{min_slope}{Minimum slope across protein}
#'   \item{min_comparison_slope}{Minimum slope for comparison of non-control
#'    \emph{vs} control}
#'   \item{min_R_sq}{Minimum \eqn{R^2} across protein}
#'   \item{max_control_plateau}{Maximum curve plateau in control samples for
#'    each protein}
#'   \item{adj_pvalue}{p-value for melting point difference, FDR adjusted using
#'    \emph{Benjamini-Hochberg} correction}}
"analysis_data_4prot"

#' Quantity data for 20 Protein Melting Curves
#'
#' A dataset containing the raw proteomics data (DIA) for 20 proteins
#'
#' @format A data frame with 800 rows and 8 variables:
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
"quan_data_20prot"

#' Statistics for 20 Protein Melting Curve
#'
#' A dataset containing the statistics and parameters generated from fitting a
#' melting curve to normalised proteomics data (DIA) for the 20 proteins and
#' comparing melting curves across measurements:
#'
#' @format A data frame with 80 rows and 19 variables:
#' \describe{
#'   \item{Protein_ID}{Gene name identifying the protein}
#'   \item{Condition}{Treatment given to this sample set}
#'   \item{Replicate}{Replicate number for this sample set}
#'   \item{F_scaled}{Scaled F-statistic from NPARC analysis}
#'   \item{p_adj_NPARC}{p-value from NPARC analysis, FDR adjusted using
#'    \emph{Benjamini-Hochberg} correction}
#'   \item{`a`}{Parameter `a` from sigmoidal curve fitting}
#'   \item{`b`}{Parameter `b` from sigmoidal curve fitting}
#'   \item{melt_point}{Melting point, calculated as the model temperature
#'    resulting in 50 % of the protein destabilised}
#'   \item{infl_point}{Inflection point, calculated as the model temperature
#'    for which \eqn{f''(T_{infl}) = 0 }}
#'   \item{slope}{Sigmoidal curve slope, calculated as \eqn{f'(T_{infl}) }}
#'   \item{plateau}{The plateau of the modeled sigmoidal curve}
#'   \item{R_sq}{\eqn{R^2} for the fitted sigmoidal curve}
#'   \item{Comparison}{Measurements compared for \eqn{\Delta T_m}}
#'   \item{diff_melt_point}{Melting point difference (\eqn{\Delta T_m})}
#'   \item{min_slope}{Minimum slope across protein}
#'   \item{min_comparison_slope}{Minimum slope for comparison of non-control
#'    \emph{vs} control}
#'   \item{min_R_sq}{Minimum \eqn{R^2} across protein}
#'   \item{max_control_plateau}{Maximum curve plateau in control samples for
#'    each protein}
#'   \item{adj_pvalue}{p-value for melting point difference, FDR adjusted using
#'    \emph{Benjamini-Hochberg} correction}}
"analysis_data_20prot"


