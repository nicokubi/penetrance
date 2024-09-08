
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Penetrance R package

An R package for the estimation of age-specific penetrance for complex
family-based studies in a format compatible with with the PanelPRO R
package.

## Motivation

Accurate estimation of age-specific penetrance is essential for
assessing cancer risk in individuals with pathogenic genetic variants
(PGVs). Penetrance refers to the probability that an individual carrying
a specific genetic variant will develop the associated trait, such as
cancer. Estimating this probability is a crucial step in clinical
decision-making and personalized risk assessment for hereditary (cancer)
syndromes.

The package leverages Mendelian inheritance models, which are widely
used in family-based genetic studies to assess how genetic variants are
passed down through generations. These models typically involve a —an
individual for whom family history and genetic data are collected. The
proband serves as the starting point for mapping out the family’s
genetic structure, including relationships and phenotypic traits, such
as cancer diagnoses. Family data, including cancer occurrence, ages of
diagnosis, and genetic test results, are collected for the proband and
their relatives. Using this data, Mendelian models compute the
likelihood of certain genetic configurations and disease outcomes based
on inheritance patterns.

The core methodology in the package relies on a four-parameter Weibull
distribution to model age-specific penetrance. Estimation is performed
using a Bayesian framework with Markov Chain Monte Carlo (MCMC) methods,
allowing the package to provide robust and flexible penetrance
estimates. Through this approach, the package models the likelihood of
cancer occurrence across family members, even when some genotypic
information is missing or incomplete, which is common in real-world
studies.

The package also incorporates prior knowledge into the estimation
process, enabling users to specify default, custom, or study-based prior
distributions. By employing the Elston-Stewart peeling algorithm, the
package efficiently calculates likelihoods across family pedigrees,
ensuring scalability and accuracy, even in large datasets.

By providing user-friendly functions for data input, prior
specification, and estimation, the package equips researchers and
clinicians with a powerful tool for estimating cancer risk in complex
family-based studies. This empowers informed decision-making and
preventive strategies in hereditary cancer syndromes, where
understanding the genetic basis of risk is critical for patient care.

## Installation

To install, use

    git clone git@github.com:nicokubi/penetrance.git

Open the source directory as new R project and install the package with

    devtools::install()

or directly in R studio

    devtools::install_github("https://github.com/nicokubi/penetrance")

## Quick-start guide

This following is a quick-start guide for basic usage of the package.
For greater detail on options, please refer to the other articles.

The primary function in the package is `PenEstim`. The package workflow
includes three main parts: user input, including family data in the form
of pedigrees and specification for the penetrance estimation, the
estimation of the posterior distribution using the MCMC approach, and
the outputting of the results in the form of the samples from the
approximated posterior distribution, i.e. the estimated penetrance
function.

``` r
library(penetrance)
```

### Pedigree

The user must specify the `pedigree` argument as a data frame with the
following columns:

-   `ID`: A numeric value representing the unique identifier for each
    individual. There should be no duplicated entries.
-   `Sex`: A numeric value where `0` indicates female and `1` indicates
    male. Missing entries are not currently supported.
-   `MotherID`: A numeric value representing the unique identifier for
    an individual’s mother.
-   `FatherID`: A numeric value representing the unique identifier for
    an individual’s father.
-   `isProband`: A numeric value where `1` indicates the individual is a
    proband and `0` otherwise.
-   `CurAge`: A numeric value indicating the age of censoring (current
    age if the person is alive or age at death if the person is
    deceased). Allowed ages range from `1` to `94`.
-   `isAff`: A numeric value indicating the affection status of cancer,
    with `1` for diagnosed individuals and `0` otherwise. Missing
    entries are not supported.
-   `Age`: A numeric value indicating the age of cancer diagnosis,
    encoded as `NA` if the individual was not diagnosed. Allowed ages
    range from `1` to `94`.
-   `Geno`: A column for germline testing or tumor marker testing
    results. Positive results should be coded as `1`, negative results
    as `0`, and unknown results as `NA` or left empty.

### Model specification

There are a few ways in which a user can specify how the estimation
approach is run. Available options are:

``` r
#' @param pedigree A data frame containing the pedigree data in the required format. It should include the following columns:
#'   - `PedigreeID`:  A numeric value representing the unique identifier for each family There should be no duplicated entries.
#'   - `ID`: A numeric value representing the unique identifier for each individual. There should be no duplicated entries.
#'   - `Sex`: A numeric value where `0` indicates female and `1` indicates male. Missing entries are not currently supported.
#'   - `MotherID`: A numeric value representing the unique identifier for an individual's mother.
#'   - `FatherID`: A numeric value representing the unique identifier for an individual's father.
#'   - `isProband`: A numeric value where `1` indicates the individual is a proband and `0` otherwise.
#'   - `CurAge`: A numeric value indicating the age of censoring (current age if the person is alive or age at death if the person is deceased). Allowed ages range from `1` to `94`.
#'   - `isAff`: A numeric value indicating the affection status of cancer, with `1` for diagnosed individuals and `0` otherwise. Missing entries are not supported.
#'   - `Age`: A numeric value indicating the age of cancer diagnosis, encoded as `NA` if the individual was not diagnosed. Allowed ages range from `1` to `94`.
#'   - `Geno`: A column for germline testing or tumor marker testing results. Positive results should be coded as `1`, negative results as `0`, and unknown results as `NA` or left empty.
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
#' @param plot_trace Boolean indicating whether to include trace plots in the output. Default is TRUE.
#' @param penetrance_plot Boolean indicating whether to include penetrance plots in the output. Default is TRUE.
#' @param penetrance_plot_pdf Boolean indicating whether to include PDF plots in the output. Default is TRUE.
#' @param probCI Probability level for confidence intervals in penetrance plots. Default is 0.95.
```

### Prior Specification

Penetrance provides the user with a flexible approach to prior
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
