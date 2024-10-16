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
    if (sex_specific) {
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
    is_male <- (sex[idx] == 1)
    rel_prob <- relationship_probs[idx]

    if (!is.na(rel_prob) && runif(1) < rel_prob) {
      # Weibull distribution for affected individuals
      imputed_ages[idx] <- get_valid_age(function() {
        if (sex_specific) {
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
        if (sex_specific) {
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

#' Calculate Degree of Relationship in Pedigree
#'
#' This function calculates the degree of relationship for individuals in a dataset based on their pedigree.
#'
#' @param data A data frame containing the data with columns 'individual', 'father', 'mother', 'sex', 'aff', 'family', and 'isProband'.
#' @return The data frame with an additional column 'degree_of_relationship' indicating the degree of relationship for each individual.
#' 
calcPedDegree <- function(data) {
  # Create a copy of the data to avoid modifying the original data directly
  data_copy <- data
  
  # Replace 0 with NA in the mother and father columns in the copy
  data_copy$mother[data_copy$mother == 0] <- NA
  data_copy$father[data_copy$father == 0] <- NA
  
  # Ensure both parents are NA if one is missing
  data_copy$father[is.na(data_copy$father) != is.na(data_copy$mother)] <- NA
  data_copy$mother[is.na(data_copy$father) != is.na(data_copy$mother)] <- NA
  
  # Create the pedigree object
  ped <- pedigree(
    id = data_copy$individual,
    dadid = data_copy$father,
    momid = data_copy$mother,
    sex = data_copy$sex,
    affected = data_copy$aff,
    famid = data_copy$family
  )
  
  # Calculate the kinship matrix
  kin_matrix <- kinship(ped)
  
  # Function to calculate the degree of relationship
  calculate_degree <- function(proband_id, family_kin_matrix) {
    degrees <- rep(NA, nrow(family_kin_matrix))
    names(degrees) <- rownames(family_kin_matrix)
    for (i in 1:nrow(family_kin_matrix)) {
      if (i == proband_id) {
        degrees[i] <- 0
      } else {
        kin_value <- family_kin_matrix[proband_id, i]
        degrees[i] <- kin_value * 2
      }
    }
    return(degrees)
  }
  
  # Initialize a column for degrees of relationship in the copy
  data_copy$degree_of_relationship <- NA
  
  # Iterate through each family to calculate the degree of relationship
  families <- unique(data_copy$family)
  
  for (fam in families) {
    family_data <- data_copy[data_copy$family == fam, ]
    proband_actual_id <- family_data$individual[family_data$isProband == 1 & !is.na(family_data$isProband)]
    
    if (length(proband_actual_id) == 1) {
      # Subset the kinship matrix for the current family
      family_ids <- family_data$individual
      family_kin_matrix <- kin_matrix[family_ids, family_ids]
      
      if (!is.null(nrow(family_kin_matrix))) {
        # Find the row index of the proband in the family kinship matrix
        proband_index <- which(family_ids == proband_actual_id)
        
        # Calculate degrees of relationship for this family
        degrees_of_relationship <- calculate_degree(proband_index, family_kin_matrix)
        
        # Update the copy of the main data frame
        data_copy$degree_of_relationship[data_copy$family == fam] <- degrees_of_relationship
      } else {
        # If the kinship matrix is invalid, set the degree_of_relationship to NA
        data_copy$degree_of_relationship[data_copy$family == fam] <- NA
      }
    } else {
      # If no unique proband is found, set the degree_of_relationship to NA
      data_copy$degree_of_relationship[data_copy$family == fam] <- NA
    }
  }
  
  # Add the new column to the original data
  data$degree_of_relationship <- data_copy$degree_of_relationship
  
  return(data)
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
  if (sex == 1) { # Male
    density_data <- if(tested) empirical_density$male_tested else empirical_density$male_untested
  } else if (sex == 2) { # Female
    density_data <- if(tested) empirical_density$female_tested else empirical_density$female_untested
  } else {
    stop("Invalid sex value.")
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
    is_male <- (sex[idx] == 1)
    is_tested <- tested[idx]

    # Draw age from the appropriate empirical distribution
    imputed_age <- round(drawEmpirical(empirical_density, ifelse(is_male, 1, 2), is_tested))

    # Ensure the imputed age is within valid range
    imputed_ages[idx] <- max(1, min(imputed_age, max_age))
  }

  # Assign imputed ages back to the data
  data$age[na_indices] <- imputed_ages

  return(data)
}