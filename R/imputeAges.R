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
#' @param empirical_density A density object or list of density objects containing the empirical density of ages, 
#'                          possibly sex-specific. If sex-specific, it should be a list with elements 'male' and 'female'.
#' @param max_age Integer, the maximum age considered in the analysis.
#' @param sex_specific Logical, indicating whether the imputation should be sex-specific. Default is TRUE.
#'
#' @return The data frame with imputed ages.
#'
imputeAges <- function(data, na_indices, baseline_male = NULL, baseline_female = NULL, 
                       alpha_male = NULL, beta_male = NULL, delta_male = NULL,
                       alpha_female = NULL, beta_female = NULL, delta_female = NULL, 
                       baseline = NULL, alpha = NULL, beta = NULL, delta = NULL,
                       empirical_density, max_age, sex_specific = TRUE) {
  
  for (i in na_indices) {
    # Ensure 'i' is within the valid range
    if (i > nrow(data) || i < 1) {
      stop(paste("Index", i, "is out of bounds for the data frame."))
    }
    
    valid_age <- FALSE
    while (!valid_age) {
      u <- runif(1)
      relationship_prob <- as.numeric(data$degree_of_relationship[i])
      
      if (is.na(relationship_prob)) {
        warning(paste("Invalid degree of relationship for individual at index", i, ". Using empirical distribution."))
        # Use drawEmpirical with sex if sex-specific
        if (sex_specific) {
          data$age[i] <- round(drawEmpirical(empirical_density, data$sex[i]))
        } else {
          data$age[i] <- round(drawEmpirical(empirical_density))
        }
      } else {
        if (data$aff[i] == 1) {
          if (runif(1) < relationship_prob) {
            if (sex_specific) {
              if (data$sex[i] == 1) { # Male
                data$age[i] <- round(delta_male + beta_male * (-log(1 - u))^(1 / alpha_male))
              } else if (data$sex[i] == 2) { # Female
                data$age[i] <- round(delta_female + beta_female * (-log(1 - u))^(1 / alpha_female))
              }
            } else {
              # Non-sex-specific imputation
              data$age[i] <- round(delta + beta * (-log(1 - u))^(1 / alpha))
            }
          } else {
            if (sex_specific) {
              data$age[i] <- round(ifelse(data$sex[i] == 1, drawBaseline(baseline_male), drawBaseline(baseline_female)))
            } else {
              # Non-sex-specific baseline imputation
              data$age[i] <- round(drawBaseline(baseline)) # Use baseline directly
            }
          }
        } else {
          # Use drawEmpirical with sex if sex-specific
          if (sex_specific) {
            data$age[i] <- round(drawEmpirical(empirical_density, data$sex[i]))
          } else {
            data$age[i] <- round(drawEmpirical(empirical_density))
          }
        }
      }
      
      if (!is.na(data$age[i]) && data$age[i] >= 1 && data$age[i] <= max_age) {
        valid_age <- TRUE
      }
    }
  }
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
    proband_actual_id <- family_data$individual[family_data$isProband == 1]
    
    if (length(proband_actual_id) == 1) {
      # Subset the kinship matrix for the current family
      family_ids <- family_data$individual
      family_kin_matrix <- kin_matrix[family_ids, family_ids]
      
      if (is.matrix(family_kin_matrix) && !is.null(nrow(family_kin_matrix))) {
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
#' This function calculates the empirical density for ages of non-affected individuals in a dataset.
#'
#' @param data A data frame containing the data.
#' @param aff_column Character, the name of the column indicating affection status. Default is "aff".
#' @param age_column Character, the name of the column indicating ages. Default is "age".
#' @param sex_column Character, the name of the column indicating sex. Default is "sex".
#' @param n_points Integer, the number of points to use in the density estimation. Default is 10000.
#' @param sex_specific Logical, indicating whether to calculate separate densities for males and females. Default is TRUE.
#'
#' @return A density object or a list of density objects (if sex-specific) representing the empirical density of ages.
#'
calculateEmpiricalDensity <- function(data, aff_column = "aff", age_column = "age", sex_column = "sex", 
                                      n_points = 10000, sex_specific = TRUE) {
  
  if (sex_specific) {
    # Calculate empirical density separately for males and females
    empirical_density <- list(
      male = density(na.omit(subset(data, data[[aff_column]] == 0 & data[[sex_column]] == 1)[[age_column]]), n = n_points),
      female = density(na.omit(subset(data, data[[aff_column]] == 0 & data[[sex_column]] == 2)[[age_column]]), n = n_points)
    )
  } else {
    # Calculate empirical density for the entire dataset (non-sex-specific)
    non_affected_data <- subset(data, data[[aff_column]] == 0)
    cleaned_non_affected_ages <- na.omit(non_affected_data[[age_column]])
    empirical_density <- density(cleaned_non_affected_ages, n = n_points)
  }
  
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
#' This function draws ages using the inverse CDF method from empirical density data.
#'
#' @param empirical_density A density object or a list of density objects (if sex-specific) containing the empirical density of ages.
#' @param sex Numeric, the sex of the individual (1 for male, 2 for female). Required if `empirical_density` is a list.
#'
#' @return A single age value drawn from the empirical density data.
#'
drawEmpirical <- function(empirical_density, sex = NA) {
  u <- runif(1)
  
  # Use the appropriate density based on the sex input
  if (is.list(empirical_density) && !is.na(sex)) {
    if (sex == 1) { # Male
      density_data <- empirical_density$male
    } else if (sex == 2) { # Female
      density_data <- empirical_density$female
    } else {
      stop("Invalid sex value for sex-specific empirical density.")
    }
  } else {
    # Use non-sex-specific density
    density_data <- empirical_density
  }
  
  age <- approx(cumsum(density_data$y) / sum(density_data$y), density_data$x, xout = u)$y
  return(age)
}