
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PenEstim R package

An R package for the estimation of age-specific penetrance for complex
family-based studies in a format compatible with with the PanelPRO R
package.

## Installation

To install, use

    git clone git@github.com:bayesmendel/PenEstim.git

Open the source directory as new R project and install the package with

    devtools::install()

## Quick-start guide

This following is a quick-start guide for basic usage of the package.
For greater detail on options, please refer to the other vignettes
[using_PenEstim]() and documentation.

The primary function in the package is `PenEstim`. The package workflow
includes three main parts: user input, including family data in the form
of pedigrees and specification for the penetrance estimation, the
estimation of the posterior distribution using the MCMC approach, and
the outputting of the results in the form of the samples from the
approximated posterior distribution, i.e. the estimated penetrance
function.

``` r
library(PenEstim)
```

### Pedigree

The user must specify the `pedigree` argument as a data frame with the
following columns: - `ID`: A numeric value representing the unique
identifier for each individual. There should be no duplicated entries. -
`Sex`: A numeric value where `0` indicates female and `1` indicates
male. Missing entries are not currently supported. - `MotherID`: A
numeric value representing the unique identifier for an individual’s
mother. - `FatherID`: A numeric value representing the unique identifier
for an individual’s father. - `isProband`: A numeric value where `1`
indicates the individual is a proband and `0` otherwise. - `CurAge`: A
numeric value indicating the age of censoring (current age if the person
is alive or age at death if the person is deceased). Allowed ages range
from `1` to `94`. - `isAff`: A numeric value indicating the affection
status of cancer, with `1` for diagnosed individuals and `0` otherwise.
Missing entries are not supported. - `Age`: A numeric value indicating
the age of cancer diagnosis, encoded as `NA` if the individual was not
diagnosed. Allowed ages range from `1` to `94`. - `Geno`: A column for
germline testing or tumor marker testing results. Positive results
should be coded as `1`, negative results as `0`, and unknown results as
`NA` or left empty.

### Model specification

There are a few ways in which a user can specify how the estimation
approach is run. Available options are:

``` r
#' @param twins Identical twins or triplets in the family can be specifed. For example, to indicate that `ora024` and `ora027` are identical twins, and so are `aey063` and `aey064`, 
#' then we can use the following as the twins arguement: twins <- list(c("ora024", "ora027"), c("aey063", "aey064"))
#' @param n_chains Number of chains for parallel computation. Default is 1.
#' @param n_iter_per_chain Number of iterations for each chain. Default is 10000.
#' @param ncores Number of cores for parallel computation. Default is 6.
#' @param baseline_data Data for the baseline risk estimates (probability of developing cancer), such as population-level risk from a cancer registry. Default is the allele frequency for MLH1 from the PanelPRO database.
#' @param max_age Maximum age considered for analysis. Default is 94.
#' @param removeProband Logical, indicating whether to remove probands from the analysis. Default is FALSE.
#' @param ageImputation Logical, indicating whether to perform age imputation. Default is FALSE.
#' @param median_max Boolean indicating whether to use the baseline median age or max_age as an upper bound for the median proposal. Default is TRUE.
#' @param BaselineNC Boolean indicating that the non-carrier penetrance is assumed to be the baseline penetrance. Default is TRUE.
#' @param var Vector of variances for the proposal distribution. Default is c(0.1, 0.1, 2, 2, 5, 5, 5, 5).
#' @param burn_in Fraction of results to discard as burn-in (0 to 1). Default is 0 (no burn-in).
#' @param thinning_factor Factor by which to thin the results. Default is 1 (no thinning).
#' @param distribution_data Data for generating prior distributions.
#' @param af Allele frequency for the risk allele. Default is 0.0001.
#' @param max_penetrance Maximum penetrance considered for analysis. Default is 1.
#' @param sample_size Optional sample size for distribution generation.
#' @param ratio Optional ratio parameter for distribution generation.
#' @param prior_params Parameters for prior distributions.
#' @param risk_proportion Proportion of risk for distribution generation.
#' @param summary_stats Boolean indicating whether to include summary statistics in the output. Default is TRUE.
#' @param rejection_rates Boolean indicating whether to include rejection rates in the output. Default is TRUE.
#' @param density_plots Boolean indicating whether to include density plots in the output. Default is TRUE.
#' @param penetrance_plot Boolean indicating whether to include penetrance plots in the output. Default is TRUE.
#' @param probCI Probability level for confidence intervals in penetrance plots. Default is 0.95.
```

### Prior Specification

PenEstim provides the user with a flexible approach to prior
specification, balancing customization with an easy-to-use workflow. In
addition to providing default prior distributions, the package allows
users to customize the priors by including existing penetrance estimates
or prior knowledge. The following settings for the prior distribution
specification are available:
### Additional User Inputs

-   The `PenEstim` function takes baseline age-specific probabilitie of
    developing cancer as as input `baseline_data`. In the default
    setting with `BaselineNC = TRUE` this baseline is assumed to reflect
    the non-carrier penetrance. For rare mutations this is considered a
    reasonable assumption. The baseline_data must be a data frame with
    baseline data for females and males.

-   Furthermore the specification of the allele frequency of the
    pathogenic variant (`af`) is required.

-   The PenEsim also includes an option for automatic age imputation
    `AgeImputation`. We apply an age imputation as part of the MCMC
    routine. The imputation of ages is performed based on the
    individual’s affected status ($aff$), sex ($sex$), and their degree
    of relationship to the proband who is a carrier of the PV. For
    greater detail on the age imputation approach see documentation.

-   For the likelihood calculation monozygous twins can be specified
    using the `twins` arguement.

``` r
twins <- list(c("ora024", "ora027"), c("aey063", "aey064"))
```
