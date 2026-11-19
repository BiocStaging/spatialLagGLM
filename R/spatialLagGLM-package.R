#' spatialLagGLM: Spatial Lag Regression for Spatial Transcriptomics
#'
#' Fits, for every gene, a quasi-Poisson regression with a library-size
#' offset, user-defined design terms and a spatially lagged covariate, and
#' stores the results in the \code{rowData} of a
#' \code{\link[SpatialExperiment]{SpatialExperiment}}.
#'
#' @keywords internal
#' @useDynLib spatialLagGLM, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom methods is
#' @importFrom utils globalVariables
"_PACKAGE"
