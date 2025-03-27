#' Processed Family Data
#'
#' A dataset containing processed information about the first simulated 130 families.
#' These families are referenced in the vigniette simulation_study.Rmd and using_penetrance.Rmd
#' The user must specify the `pedigree` argument as a data frame which contains the
#' family data (see `test_fam`). The family data must be in the correct format with the following columns:
#' 
#' \describe{
#'   \item{\code{ID}}{A numeric value representing the unique identifier for each individual. There should be no duplicated entries.}
#'   \item{\code{Sex}}{A numeric value where \code{0} indicates female and \code{1} indicates male. Missing entries are not currently supported.}
#'   \item{\code{MotherID}}{A numeric value representing the unique identifier for an individual's mother.}
#'   \item{\code{FatherID}}{A numeric value representing the unique identifier for an individual's father.}
#'   \item{\code{isProband}}{A numeric value where \code{1} indicates the individual is a proband and \code{0} otherwise.}
#'   \item{\code{CurAge}}{A numeric value indicating the age of censoring (current age if the person is alive or age at death if the person is deceased). Allowed ages range from \code{1} to \code{94}.}
#'   \item{\code{isAff}}{A numeric value indicating the affection status of cancer, with \code{1} for diagnosed individuals and \code{0} otherwise. Missing entries are not supported.}
#'   \item{\code{Age}}{A numeric value indicating the age of cancer diagnosis, encoded as \code{NA} if the individual was not diagnosed. Allowed ages range from \code{1} to \code{94}.}
#'   \item{\code{Geno}}{A column for germline testing or tumor marker testing results. Positive results should be coded as \code{1}, negative results as \code{0}, and unknown results as \code{NA} or left empty.}
#' }
#' 
#' @format A list of processed family data.
#' @source Generated for package example
"simulated_families"

#' Processed Family Data
#'
#' A dataset containing processed information about the first simulated 130 families.
#' These families are referenced in the vigniette simulation_study_real.Rmd
#' The user must specify the `pedigree` argument as a data frame which contains the
#' family data (see `test_fam`). The family data must be in the correct format with the following columns:
#' 
#' \describe{
#'   \item{\code{ID}}{A numeric value representing the unique identifier for each individual. There should be no duplicated entries.}
#'   \item{\code{Sex}}{A numeric value where \code{0} indicates female and \code{1} indicates male. Missing entries are not currently supported.}
#'   \item{\code{MotherID}}{A numeric value representing the unique identifier for an individual's mother.}
#'   \item{\code{FatherID}}{A numeric value representing the unique identifier for an individual's father.}
#'   \item{\code{isProband}}{A numeric value where \code{1} indicates the individual is a proband and \code{0} otherwise.}
#'   \item{\code{CurAge}}{A numeric value indicating the age of censoring (current age if the person is alive or age at death if the person is deceased). Allowed ages range from \code{1} to \code{94}.}
#'   \item{\code{isAff}}{A numeric value indicating the affection status of cancer, with \code{1} for diagnosed individuals and \code{0} otherwise. Missing entries are not supported.}
#'   \item{\code{Age}}{A numeric value indicating the age of cancer diagnosis, encoded as \code{NA} if the individual was not diagnosed. Allowed ages range from \code{1} to \code{94}.}
#'   \item{\code{Geno}}{A column for germline testing or tumor marker testing results. Positive results should be coded as \code{1}, negative results as \code{0}, and unknown results as \code{NA} or left empty.}
#' }
#' 
#' @format A list of processed family data.
#' @source Generated for package example
"test_fam2"

#' Simulated Output Data
#'
#' This dataset contains the simulated output data for the penetrance package.
#'
#' @docType data
#' @name out_sim
#' @usage data(out_sim)
#' @format A list with the following components:
#' \describe{
#'   \item{summary_stats}{A data frame with 18000 observations of 8 variables:
#'     \describe{
#'       \item{Median_Male}{numeric, Median value for males}
#'       \item{Median_Female}{numeric, Median value for females}
#'       \item{Threshold_Male}{numeric, Threshold value for males}
#'       \item{Threshold_Female}{numeric, Threshold value for females}
#'       \item{First_Quartile_Male}{numeric, First quartile value for males}
#'       \item{First_Quartile_Female}{numeric, First quartile value for females}
#'       \item{Asymptote_Male}{numeric, Asymptote value for males}
#'       \item{Asymptote_Female}{numeric, Asymptote value for females}
#'     }
#'   }
#'   \item{density_plots}{A list of 1 element, mfrow: integer vector of length 2}
#'   \item{trace_plots}{A list of 1 element, mfrow: integer vector of length 2}
#'   \item{penetrance_plot}{A list of 2 elements: rect and text}
#'   \item{pdf_plots}{A list of 2 elements: rect and text}
#'   \item{combined_chains}{A list of 19 numeric vectors with 18000 elements each}
#'   \item{results}{A list of 1 element which is a list of 24 elements, each with 18000 elements}
#'   \item{data}{A data frame with 4727 observations of 9 variables:
#'     \describe{
#'       \item{individual}{integer, Individual ID}
#'       \item{isProband}{numeric, Indicator if the individual is a proband}
#'       \item{family}{integer, Family ID}
#'       \item{mother}{numeric, Mother's ID}
#'       \item{father}{numeric, Father's ID}
#'       \item{aff}{numeric, Affected status}
#'       \item{sex}{numeric, Sex of the individual}
#'       \item{age}{numeric, Age of the individual}
#'       \item{geno}{character, Genotype}
#'     }
#'   }
#' }
#' @examples
#' data(out_sim)
#' head(out_sim$summary_stats)
"out_sim"

#' Default Baseline Data
#'
#' This dataset contains age-specific cancer penetrance rates for both females and males.
#' As an example, the data is derived from the SEER program for colorectal cancer females and males. 
#'
#' @format A data frame with 94 rows and 3 variables:
#' \describe{
#'   \item{Age}{Age in years (1 to 94)}
#'   \item{Female}{Penetrance rate for females}
#'   \item{Male}{Penetrance rate for males}
#' }
#' @export
baseline_data_default <- data.frame(
  Age = 1:94,
  Female = c(
    0.00000005, 0.00000009, 0.00000021, 0.00000041, 0.00000061, 0.00000081, 0.00000113,
    0.00000225, 0.00000350, 0.00000474, 0.00000599, 0.00000728, 0.00000886, 0.00001049,
    0.00001212, 0.00001375, 0.00001531, 0.00001642, 0.00001746, 0.00001849, 0.00001953,
    0.00002073, 0.00002289, 0.00002521, 0.00002753, 0.00002985, 0.00003262, 0.00003808,
    0.00004400, 0.00004992, 0.00005583, 0.00006215, 0.00007085, 0.00007996, 0.00008906,
    0.00009817, 0.00010812, 0.00012320, 0.00013913, 0.00015506, 0.00017098, 0.00018798,
    0.00021142, 0.00023593, 0.00026044, 0.00028494, 0.00031276, 0.00036047, 0.00041150,
    0.00046251, 0.00051351, 0.00055826, 0.00056563, 0.00056677, 0.00056790, 0.00056903,
    0.00057393, 0.00060153, 0.00063290, 0.00066426, 0.00069559, 0.00072936, 0.00077779,
    0.00082863, 0.00087943, 0.00093019, 0.00098016, 0.00102557, 0.00107018, 0.00111473,
    0.00115923, 0.00120791, 0.00128200, 0.00136022, 0.00143831, 0.00151626, 0.00159429,
    0.00167343, 0.00175258, 0.00183150, 0.00191019, 0.00198761, 0.00205874, 0.00212853,
    0.00219797, 0.00226703, 0.00232237, 0.00229747, 0.00225914, 0.00222074, 0.00218226,
    0.00213527, 0.00203759, 0.00193173
  ),
  Male = c(
    0.00000004, 0.00000009, 0.00000022, 0.00000045, 0.00000068, 0.00000091, 0.00000118,
    0.00000174, 0.00000235, 0.00000296, 0.00000356, 0.00000424, 0.00000539, 0.00000661,
    0.00000783, 0.00000905, 0.00001025, 0.00001134, 0.00001240, 0.00001346, 0.00001452,
    0.00001587, 0.00001894, 0.00002231, 0.00002567, 0.00002904, 0.00003262, 0.00003753,
    0.00004265, 0.00004778, 0.00005290, 0.00005860, 0.00006771, 0.00007738, 0.00008706,
    0.00009673, 0.00010758, 0.00012550, 0.00014460, 0.00016369, 0.00018278, 0.00020326,
    0.00023204, 0.00026221, 0.00029237, 0.00032253, 0.00035783, 0.00042408, 0.00049547,
    0.00056684, 0.00063818, 0.00070320, 0.00073034, 0.00075116, 0.00077195, 0.00079273,
    0.00081720, 0.00086389, 0.00091425, 0.00096456, 0.00101482, 0.00106685, 0.00112977,
    0.00119444, 0.00125904, 0.00132355, 0.00138738, 0.00144751, 0.00150694, 0.00156625,
    0.00162545, 0.00168775, 0.00176926, 0.00185381, 0.00193813, 0.00202223, 0.00210529,
    0.00218325, 0.00226009, 0.00233659, 0.00241273, 0.00248511, 0.00253680, 0.00258475,
    0.00263232, 0.00267951, 0.00271651, 0.00269447, 0.00266248, 0.00263032, 0.00259802,
    0.00255423, 0.00244243, 0.00231964
  )
)