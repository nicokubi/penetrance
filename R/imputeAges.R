#' Impute Ages Based on Affection Status and Sex
#'
#' This function imputes ages for individuals in a dataset based on their affection
#' status and sex using either Weibull, baseline, or empirical distribution.
#'
#' @param data A data frame containing the individual data, including columns for age, sex, and affection status.
#' @param na_indices A vector of indices indicating the rows in the data where ages need to be imputed.
#' @param baseline_male A data frame containing baseline data for males, with columns 'cum_prob' and 'age'.
#' @param baseline_female A data frame containing baseline data for females, with columns 'cum_prob' and 'age'.
#' @param alpha_male Numeric, shape parameter for the Weibull distribution for males.
#' @param beta_male Numeric, scale parameter for the Weibull distribution for males.
#' @param delta_male Numeric, location parameter for the Weibull distribution for males.
#' @param alpha_female Numeric, shape parameter for the Weibull distribution for females.
#' @param beta_female Numeric, scale parameter for the Weibull distribution for females.
#' @param delta_female Numeric, location parameter for the Weibull distribution for females.
#' @param baseline A data frame containing baseline data (used for non-sex-specific analysis) with columns 'cum_prob' and 'age'.
#' @param alpha Numeric, shape parameter for the Weibull distribution (used for non-sex-specific analysis).
#' @param beta Numeric, scale parameter for the Weibull distribution (used for non-sex-specific analysis).
#' @param delta Numeric, location parameter for the Weibull distribution (used for non-sex-specific analysis).
#' @param max_age Integer, the maximum age considered in the analysis.
#' @param sex_specific Logical, indicating whether the imputation should be sex-specific. Default is TRUE.
#' @param max_attempts Integer, the maximum number of attempts to get a valid age. Default is 100.
#'
#' @return The data frame with imputed ages.
#'
imputeAges <- function(data, na_indices, baseline_male = NULL, baseline_female = NULL,
                       alpha_male = NULL, beta_male = NULL, delta_male = NULL,
                       alpha_female = NULL, beta_female = NULL, delta_female = NULL,
                       baseline = NULL, alpha = NULL, beta = NULL, delta = NULL,
                       max_age, sex_specific = TRUE, max_attempts = 100,
                       geno_freq, trans, lik)  {
  # Change the data frame to include the 'indiv' column and the 'mother' and 'father' columns
  data$indiv <- paste(data$family, data$individual, sep = "-")
  data$mother <- paste(data$family, data$mother, sep = "-")
  data$father <- paste(data$family, data$father, sep = "-") 
  # Remove the original 'individual', 'mother', and 'father' columns
  data$individual <- NULL

  # Move the 'indiv' column to the beginning of the dataframe
  data <- data[, c("indiv", setdiff(names(data), "indiv"))]
  
  # Extract necessary columns for faster access
  aff <- data$aff[na_indices]
  sex <- data$sex[na_indices]

  # Calculate median ages for fallback option
  median_ages <- list(
    male_affected = median(data$age[data$sex == 1 & data$aff == 1], na.rm = TRUE),
    female_affected = median(data$age[data$sex == 2 & data$aff == 1], na.rm = TRUE),
    all_affected = median(data$age[data$aff == 1], na.rm = TRUE)
  )

  # Function to get a valid age
  get_valid_age <- function(age_func, is_male) {
    for (attempt in 1:max_attempts) {
      age <- age_func()
      if (!is.na(age) && age >= 1 && age <= max_age) {
        return(age)
      }
    }
    # Fallback to median age if max attempts reached
    if (is.na(is_male)) {
      return(median_ages$all_affected)
    } else if (sex_specific) {
      return(if (is_male) median_ages$male_affected else median_ages$female_affected)
    } else {
      return(median_ages$all_affected)
    }
  }

  # Create a vector to store the imputed ages
  imputed_ages <- numeric(length(na_indices))

  # Calculate genotype probabilities for all individuals using the clipp function
  all_genotype_probs <- lapply(data$indiv, function(target) {
    genotype_probabilities(target, data, geno_freq, trans, lik)
  })

  # Extract the second column (relationship probabilities) for each individual
  relationship_probs <- sapply(all_genotype_probs, function(x) x[2])

    # Select only the relationship probabilities for individuals with NA ages
    relationship_probs <- relationship_probs[na_indices]

  for (idx in seq_along(na_indices)) {
    is_male <- if (is.na(sex[idx])) NA else (sex[idx] == 1)
    rel_prob <- relationship_probs[idx]

    if (!is.na(rel_prob) && runif(1) < rel_prob) {
      # Weibull distribution for affected individuals
      imputed_ages[idx] <- get_valid_age(function() {
        if (is.na(is_male)) {
          # Randomly choose male or female parameters
          params <- sample(
            list(
              list(delta = delta_male, beta = beta_male, alpha = alpha_male),
              list(delta = delta_female, beta = beta_female, alpha = alpha_female)
            ), 1)[[1]]
          round(params$delta + params$beta * (-log(1 - runif(1)))^(1 / params$alpha))
        } else if (sex_specific) {
          if (is_male) {
            round(delta_male + beta_male * (-log(1 - runif(1)))^(1 / alpha_male))
          } else {
            round(delta_female + beta_female * (-log(1 - runif(1)))^(1 / alpha_female))
          }
        } else {
          round(delta + beta * (-log(1 - runif(1)))^(1 / alpha))
        }
      }, is_male)
    } else {
      # Baseline for affected if not using Weibull
      imputed_ages[idx] <- get_valid_age(function() {
        if (is.na(is_male)) {
          # Randomly choose male or female baseline
          chosen_baseline <- sample(list(baseline_male, baseline_female), 1)[[1]]
          round(drawBaseline(chosen_baseline))
        } else if (sex_specific) {
          if (is_male) round(drawBaseline(baseline_male)) else round(drawBaseline(baseline_female))
        } else {
          round(drawBaseline(baseline))
        }
      }, is_male)
    }
  }

  # Assign imputed ages back to the data
  data$age[na_indices] <- imputed_ages

  # Reformat individual, mother, and father columns back to original format
  data$individual <- sapply(strsplit(data$indiv, "-"), `[`, 2)
  data$mother <- sapply(strsplit(data$mother, "-"), `[`, 2)
  data$father <- sapply(strsplit(data$father, "-"), `[`, 2)
  
  # Remove the 'indiv' column
  data$indiv <- NULL
  
  # Reorder columns to match original format
  original_order <- c("family", "individual", "father", "mother", "sex", "aff", "age", "geno", "isProband", "degree_of_relationship")
  data <- data[, original_order]

  return(data)
}

#' Initialize Ages Using a Uniform Distribution
#'
#' This function initializes ages for individuals in a dataset using a uniform distribution.
#'
#' @param data A data frame containing the data.
#' @param threshold The minimum age value for the uniform distribution.
#' @param max_age The maximum age value for the uniform distribution.
#' @return A list containing the updated data frame and the indices of the rows where ages were initialized.
#' 
imputeAgesInit <- function(data, threshold, max_age) {
  na_indices <- which(is.na(data$age))
  data$age[na_indices] <- runif(length(na_indices), threshold, max_age)
  return(list(data = data, na_indices = na_indices))
}

#' Calculate Empirical Density for Non-Affected Individuals
#'
#' This function calculates the empirical density for ages of non-affected individuals in a dataset,
#' differentiating by sex and whether the individual was tested (has a non-NA 'geno' value).
#'
#' @param data A data frame containing the data.
#' @param aff_column Character, the name of the column indicating affection status. Default is "aff".
#' @param age_column Character, the name of the column indicating ages. Default is "age".
#' @param sex_column Character, the name of the column indicating sex. Default is "sex".
#' @param geno_column Character, the name of the column indicating genotype. Default is "geno".
#' @param n_points Integer, the number of points to use in the density estimation. Default is 10000.
#'
#' @return A list of density objects representing the empirical density of ages for different groups.
#'
calculateEmpiricalDensity <- function(data, aff_column = "aff", age_column = "age", sex_column = "sex", 
                                      geno_column = "geno", n_points = 10000) {
  
  # Helper function to check if geno is missing or empty
  is_geno_missing <- function(x) is.na(x) | x == ""
  
  empirical_density <- list(
    male_tested = density(na.omit(subset(data, data[[aff_column]] == 0 & data[[sex_column]] == 1 & !is_geno_missing(data[[geno_column]]))[[age_column]]), n = n_points),
    male_untested = density(na.omit(subset(data, data[[aff_column]] == 0 & data[[sex_column]] == 1 & is_geno_missing(data[[geno_column]]))[[age_column]]), n = n_points),
    female_tested = density(na.omit(subset(data, data[[aff_column]] == 0 & data[[sex_column]] == 2 & !is_geno_missing(data[[geno_column]]))[[age_column]]), n = n_points),
    female_untested = density(na.omit(subset(data, data[[aff_column]] == 0 & data[[sex_column]] == 2 & is_geno_missing(data[[geno_column]]))[[age_column]]), n = n_points)
  )
  
  return(empirical_density)
}


#' Draw Ages Using the Inverse CDF Method from the baseline data
#'
#' This function draws ages using the inverse CDF method from baseline data.
#'
#' @param baseline_data A data frame containing baseline data with columns 'cum_prob' and 'age'.
#' @return A single age value drawn from the baseline data.
#' 
drawBaseline <- function(baseline_data) {
  u <- runif(1)
  age <- approx(baseline_data$cum_prob, baseline_data$age, xout = u)$y
  return(age)
}

#' Draw Ages Using the Inverse CDF Method from Empirical Density
#'
#' This function draws ages using the inverse CDF method from empirical density data,
#' based on sex and whether the individual was tested.
#'
#' @param empirical_density A list of density objects containing the empirical density of ages for different groups.
#' @param sex Numeric, the sex of the individual (1 for male, 2 for female).
#' @param tested Logical, indicating whether the individual was tested (has a non-NA 'geno' value).
#'
#' @return A single age value drawn from the appropriate empirical density data.
#'
drawEmpirical <- function(empirical_density, sex, tested) {
  u <- runif(1)
  
  # Select the appropriate density based on sex and tested status
  if (is.na(sex)) {
    # Use combined male/female density or pick randomly between male/female
    density_data <- if(tested) 
        sample(list(empirical_density$male_tested, empirical_density$female_tested), 1)[[1]]
    else 
        sample(list(empirical_density$male_untested, empirical_density$female_untested), 1)[[1]]
  } else if (sex == 1) { # Male
    density_data <- if(tested) empirical_density$male_tested else empirical_density$male_untested
  } else if (sex == 2) { # Female
    density_data <- if(tested) empirical_density$female_tested else empirical_density$female_untested
  } else {
    stop("Invalid sex value: must be 1, 2, or NA")
  }
  
  age <- approx(cumsum(density_data$y) / sum(density_data$y), density_data$x, xout = u)$y
  return(age)
}
#' Impute Ages for Unaffected Individuals
#'
#' This function imputes ages for unaffected individuals in a dataset based on their sex
#' and whether they were tested, using empirical age distributions.
#'
#' @param data A data frame containing the individual data, including columns for age, sex, and geno.
#' @param na_indices A vector of indices indicating the rows in the data where ages need to be imputed.
#' @param empirical_density A list of density objects containing the empirical density of ages for different groups.
#' @param max_age Integer, the maximum age considered in the analysis.
#'
#' @return The data frame with imputed ages for unaffected individuals.
#'
imputeUnaffectedAges <- function(data, na_indices, empirical_density, max_age) {
  # Extract necessary columns for faster access
  sex <- data$sex[na_indices]
  tested <- !is.na(data$geno[na_indices])

  # Create a vector to store the imputed ages
  imputed_ages <- numeric(length(na_indices))

  for (idx in seq_along(na_indices)) {
    # Pass sex value directly (can be 1, 2, or NA)
    is_tested <- tested[idx]

    # Draw age from the appropriate empirical distribution
    imputed_age <- round(drawEmpirical(empirical_density, sex[idx], is_tested))

    # Ensure the imputed age is within valid range
    imputed_ages[idx] <- max(1, min(imputed_age, max_age))
  }

  # Assign imputed ages back to the data
  data$age[na_indices] <- imputed_ages

  return(data)
}