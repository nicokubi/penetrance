% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/outputHelpers.R
\name{plot_pdf}
\alias{plot_pdf}
\title{Plot Weibull Probability Density Function with Credible Intervals}
\usage{
plot_pdf(data, prob, max_age, sex = "NA")
}
\arguments{
\item{data}{Data frame, containing individual demographic and genetic information. Must include columns for 'sex', 'age', 'aff' (affection status), and 'geno' (genotype).}

\item{prob}{Numeric, the probability level for the credible intervals. Must be between 0 and 1.}

\item{max_age}{Integer, the maximum age considered in the analysis.}

\item{sex}{Character, specifying the sex of the individuals for the plot ("Male", "Female", or "NA" for not applicable). Default is "NA".}
}
\value{
A plot showing the Weibull PDF with credible intervals.
}
\description{
This function plots the Weibull PDF with credible intervals for the given data.
It allows for visualization of density curves for individuals based on their genetic
and demographic information.
}
