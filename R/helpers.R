#' Calculate Weibull Parameters
#'
#' This function calculates the shape (\code{alpha}) and scale (\code{beta}) parameters
#' of a Weibull distribution given the median, first quartile, and delta values.
#'
#' @param given_median The median of the data.
#' @param given_first_quartile The first quartile of the data.
#' @param delta A constant offset value.
#'
#' @return A list containing the calculated Weibull parameters, \code{alpha} and \code{beta}.
#' @export
calculate_weibull_parameters <- function(given_median, given_first_quartile, delta) {
  # Calculate alpha
  alpha <- log(-log(0.5) / -log(0.75)) / log((given_median - delta) / (given_first_quartile - delta))
  
  # Calculate beta using the median (M)
  beta <- (given_median - delta) / (-log(0.5))^(1 / alpha)
  
  return(list(alpha = alpha, beta = beta))
}

#' Validate Weibull Parameters
#'
#' This function validates the given parameters for calculating Weibull distribution.
#'
#' @param given_first_quartile The first quartile of the data.
#' @param given_median The median of the data.
#' @param threshold A constant threshold value.
#' @param asymptote A constant asymptote value (gamma).
#'
#' @return Boolean indicating whether the parameters are valid (TRUE) or not (FALSE).
#'
validate_weibull_parameters <- function(given_first_quartile, given_median, threshold, asymptote) {
  # Check for negative or zero values
  if (given_median <= 0 || given_first_quartile <= 0 || threshold < 0) {
    return(FALSE)
  }
  
  # Check if asymptote (gamma) is within the valid range (0,1)
  if (asymptote <= 0 || asymptote >= 1) {
    return(FALSE)
  }
  
  # Check if the logarithmic calculations will be valid
  if (given_first_quartile <= threshold || given_median <= threshold) {
    return(FALSE)
  }
  
  # Check if the denominator in the alpha calculation would be zero
  if ((given_first_quartile - threshold) == (given_median - threshold)) {
    return(FALSE)
  }
  
  # If all checks pass, return TRUE
  return(TRUE)
}

#' Transform Data Frame
#'
#' This function transforms a data frame from the standard format used in PanelPRO
#' into the required format which conforms to the requirements of penetrance (and clipp).
#'
#' @param df The input data frame in the usual PanelPRO format.
#'
#' @return The transformed data frame in the format required for clipp.
#' @export
transformDF <- function(df) {
  # Rename and transform columns
  df$individual <- df$ID
  df$isProband <- df$isProband
  df$family <- df$PedigreeID
  df$mother <- df$MotherID
  df$father <- df$FatherID
  df$aff <- df$isAff
  df$sex <- ifelse(df$Sex == 0, 2, df$Sex) # Convert 0s to 2s (female), keep 1s as is (male)
  
  # Apply row-wise logic to assign age based on aff column
  df$age <- apply(df, 1, function(row) {
    aff <- as.numeric(row[["isAff"]])
    if (aff == 1) {
      return(as.numeric(row[["Age"]]))
    } else {
      return(as.numeric(row[["CurAge"]]))
    }
  })
  
  # Process 'geno' column, checking for both 'geno' and 'Geno'
  if ("geno" %in% colnames(df)) {
    df$geno <- ifelse(is.na(df$geno), "", ifelse(df$geno == 1, "1/2", ifelse(df$geno == 0, "1/1", df$geno)))
  } else if ("Geno" %in% colnames(df)) {
    df$geno <- ifelse(is.na(df$Geno), "", ifelse(df$Geno == 1, "1/2", ifelse(df$Geno == 0, "1/1", df$Geno)))
  }
  
  # Select only the necessary columns
  df <- df[c("individual", "isProband", "family", "mother", "father", "aff", "sex", "age", "geno")]
  
  return(df)
}