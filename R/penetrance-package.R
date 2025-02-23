#' penetrance: A Package for Penetrance Estimation
#'
#' @description 
#' A comprehensive package for penetrance estimation in family-based studies. This package
#' implements Bayesian methods using Metropolis-Hastings algorithm for estimating age-specific
#' penetrance of genetic variants. It supports both sex-specific and non-sex-specific analyses,
#' and provides various visualization tools for examining MCMC results.
#'
#' @details
#' Key features:
#' \itemize{
#'   \item Bayesian estimation of penetrance using family-based data
#'   \item Support for sex-specific and non-sex-specific analyses
#'   \item Age imputation for missing data
#'   \item Visualization tools for MCMC diagnostics
#'   \item Integration with the clipp package for likelihood calculations
#' }
#'
#' @name penetrance
#' @aliases penetrance-package
#' @import MASS
#' @import kinship2
#' @importFrom grDevices adjustcolor dev.list dev.new rgb
#' @importFrom graphics grid hist legend lines par polygon
#' @importFrom stats acf approx cov dbeta density dunif dweibull median 
#'             na.omit pweibull qbeta quantile var
#' @importFrom clipp genotype_probabilities pedigree_loglikelihood
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL 