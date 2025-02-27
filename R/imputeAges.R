#' Impute Missing Ages in Family-Based Data
#'
#' @description
#' Imputes missing ages in family-based data using a combination of Weibull distributions
#' for affected individuals and empirical distributions for unaffected individuals.
#' The function can perform both sex-specific and non-sex-specific imputations.
#'
#' @param data A data frame containing family-based data with columns: family, individual,
#'   father, mother, sex, aff, age, geno, and isProband
#' @param na_indices Vector of indices where ages need to be imputed
#' @param baseline_male,baseline_female Data frames containing baseline age distributions for males/females
#' @param alpha_male,alpha_female Shape parameters for male/female Weibull distributions
#' @param beta_male,beta_female Scale parameters for male/female Weibull distributions
#' @param delta_male,delta_female Location parameters for male/female Weibull distributions
#' @param baseline Data frame containing overall baseline age distribution (non-sex-specific)
#' @param alpha,beta,delta Overall Weibull parameters (non-sex-specific)
#' @param max_age Maximum allowable age
#' @param sex_specific Logical; whether to use sex-specific parameters
#' @param max_attempts Maximum number of attempts for generating valid ages
#' @param geno_freq Vector of genotype frequencies
#' @param trans Transmission probabilities
#' @param lik Likelihood matrix
#'
#' @return A data frame with the following modifications:
#'   \item{age}{Updated with imputed ages for previously missing values}
#'   The rest of the data frame remains unchanged.
#' 
#' @examples
#' # Create sample data with the same structure as used in mhChain
#' data <- data.frame(
#'   family = rep(1:2, each=5),
#'   individual = rep(1:5, 2),
#'   father = c(NA,1,1,1,1, NA,6,6,6,6),
#'   mother = c(NA,2,2,2,2, NA,7,7,7,7),
#'   sex = c(1,2,1,2,1, 1,2,1,2,1),
#'   aff = c(1,0,1,0,NA, 1,0,1,0,NA),
#'   age = c(45,NA,25,NA,20, 50,NA,30,NA,22),
#'   geno = c("1/2",NA,"1/2",NA,NA, "1/2",NA,"1/2",NA,NA),
#'   isProband = c(1,0,0,0,0, 1,0,0,0,0)
#' )
#' 
#' # Initialize parameters
#' na_indices <- which(is.na(data$age))
#' geno_freq <- c(0.999, 0.001)  # Frequency of normal and risk alleles
#' trans <- matrix(c(1,0,0.5,0.5), nrow=2)  # Transmission matrix
#' lik <- matrix(1, nrow=nrow(data), ncol=2)  # Likelihood matrix
#' 
#' # Create baseline data for both sex-specific and non-sex-specific cases
#' age_range <- 20:94
#' n_ages <- length(age_range)
#' 
#' # Sex-specific baseline data
#' baseline_male <- data.frame(
#'   age = age_range,
#'   cum_prob = (1:n_ages)/n_ages * 0.8  # Male cumulative probabilities
#' )
#' 
#' baseline_female <- data.frame(
#'   age = age_range,
#'   cum_prob = (1:n_ages)/n_ages * 0.9  # Female cumulative probabilities
#' )
#' 
#' # Non-sex-specific baseline data
#' baseline <- data.frame(
#'   age = age_range,
#'   cum_prob = (1:n_ages)/n_ages * 0.85  # Overall cumulative probabilities
#' )
#' 
#' # Example with sex-specific imputation
#' imputed_data_sex <- imputeAges(
#'   data = data,
#'   na_indices = na_indices,
#'   baseline_male = baseline_male,
#'   baseline_female = baseline_female,
#'   alpha_male = 3.5,
#'   beta_male = 20,
#'   delta_male = 20,
#'   alpha_female = 3.2,
#'   beta_female = 18,
#'   delta_female = 18,
#'   max_age = 94,
#'   sex_specific = TRUE,
#'   geno_freq = geno_freq,
#'   trans = trans,
#'   lik = lik
#' )
#' 
#' # Example with non-sex-specific imputation
#' imputed_data_nosex <- imputeAges(
#'   data = data,
#'   na_indices = na_indices,
#'   baseline = baseline,
#'   alpha = 3.3,
#'   beta = 19,
#'   delta = 19,
#'   max_age = 94,
#'   sex_specific = FALSE,
#'   geno_freq = geno_freq,
#'   trans = trans,
#'   lik = lik
#' )
#' @export
imputeAges <- function(data, na_indices, baseline_male = NULL, baseline_female = NULL,
                       alpha_male = NULL, beta_male = NULL, delta_male = NULL,
                       alpha_female = NULL, beta_female = NULL, delta_female = NULL,
                       baseline = NULL, alpha = NULL, beta = NULL, delta = NULL,
                       max_age, sex_specific = TRUE, max_attempts = 100,
                       geno_freq, trans, lik) {
  
  # Store original IDs once at the start
  original_ids <- list(
    individual = data$individual,
    mother = data$mother,
    father = data$father
  )

  # Transform IDs
  data$indiv <- paste(data$family, data$individual, sep = "-")
  data$mother <- paste(data$family, data$mother, sep = "-")
  data$father <- paste(data$family, data$father, sep = "-") 
  data$individual <- NULL

  # Split indices by affection status
  affected_indices <- na_indices[data$aff[na_indices] == 1]
  unaffected_indices <- na_indices[data$aff[na_indices] == 0]

  # Calculate empirical density for unaffected individuals
  age_density <- calculateEmpiricalDensity(data, sex_specific = sex_specific)

  # Handle affected individuals
  if (length(affected_indices) > 0) {
    # Extract necessary columns for faster access
    aff <- data$aff[affected_indices]
    sex <- data$sex[affected_indices]

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
    imputed_ages <- numeric(length(affected_indices))

    # Calculate genotype probabilities for individuals with NA ages
    all_genotype_probs <- lapply(seq_along(affected_indices), function(i) {
      idx <- affected_indices[i]
      target <- data$indiv[idx]
      
      result <- tryCatch({
        # Pass the full lik matrix but use the correct index
        probs <- genotype_probabilities(target, data, geno_freq, trans, lik)
        
        if (is.null(probs) || length(probs) == 0 || all(is.na(probs))) {
          warning(sprintf("Invalid probability calculation for %s", target))
          return(c(NA, NA))
        }
        
        probs
        
      }, error = function(e) {
        warning(sprintf("Error calculating genotype probabilities for %s: %s", target, e$message))
        return(c(NA, NA))
      })
      
      return(result)
    })

    # Extract the second column (relationship probabilities)
    relationship_probs <- sapply(all_genotype_probs, function(x) x[2])

    # Handle cases where all probabilities are NA
    if (all(is.na(relationship_probs))) {
      warning("All genotype probability calculations failed, using baseline probabilities")
      relationship_probs <- rep(0, length(relationship_probs))
    }

    for (idx in seq_along(affected_indices)) {
      is_male <- if (is.na(sex[idx])) NA else (sex[idx] == 1)
      rel_prob <- relationship_probs[idx]

      if (!is.na(rel_prob) && runif(1) < rel_prob) {
        # Weibull distribution for affected individuals
        imputed_ages[idx] <- get_valid_age(function() {
          if (sex_specific) {
            if (is.na(is_male)) {
              return(NA)  # Return NA if sex is unknown in sex-specific analysis
            }
            if (is_male) {
              round(delta_male + beta_male * (-log(1 - runif(1)))^(1 / alpha_male))
            } else {
              round(delta_female + beta_female * (-log(1 - runif(1)))^(1 / alpha_female))
            }
          } else {
            # Non-sex-specific analysis: use common parameters regardless of sex
            round(delta + beta * (-log(1 - runif(1)))^(1 / alpha))
          }
        }, is_male)
      } else {
        # Baseline for affected if not using Weibull
        imputed_ages[idx] <- get_valid_age(function() {
          if (sex_specific) {
            if (is.na(is_male)) {
              return(NA)  # Return NA if sex is unknown in sex-specific analysis
            }
            if (is_male) {
              round(drawBaseline(baseline_male))
            } else {
              round(drawBaseline(baseline_female))
            }
          } else {
            # Non-sex-specific analysis: use common baseline regardless of sex
            round(drawBaseline(baseline))
          }
        }, is_male)
      }
    }

    # Assign imputed ages back to the data
    data$age[affected_indices] <- imputed_ages
  }

  # Handle unaffected individuals
  if (length(unaffected_indices) > 0) {
    sex <- data$sex[unaffected_indices]
    tested <- !is.na(data$geno[unaffected_indices])
    
    for (idx in seq_along(unaffected_indices)) {
      is_tested <- tested[idx]
      # Draw age from empirical distribution
      imputed_age <- round(drawEmpirical(age_density, sex[idx], is_tested, sex_specific = sex_specific))
      # Ensure age is within valid range
      data$age[unaffected_indices[idx]] <- max(1, min(imputed_age, max_age))
    }
  }

  # Single restoration of original format at the end
  data$individual <- original_ids$individual
  data$mother <- original_ids$mother
  data$father <- original_ids$father
  data$indiv <- NULL
  
  # Reorder columns
  original_order <- c("family", "individual", "father", "mother", "sex", "aff", "age", "geno", "isProband")
  data <- data[, original_order]

  return(data)
}

#' Initialize Age Imputation
#'
#' @description
#' Initializes the age imputation process by filling missing ages with random values
#' between a threshold and maximum age.
#'
#' @param data A data frame containing family-based data
#' @param threshold Minimum age value for initialization
#' @param max_age Maximum age value for initialization
#'
#' @return A list containing:
#'   \item{data}{The data frame with initialized ages}
#'   \item{na_indices}{Indices of missing age values}
#' @export
#' @examples
#' # Create sample data
#' data <- data.frame(
#'   family = c(1, 1),
#'   individual = c(1, 2),
#'   father = c(NA, 1),
#'   mother = c(NA, NA),
#'   sex = c(1, 2),
#'   aff = c(1, 0),
#'   age = c(NA, NA),
#'   geno = c("1/2", NA),
#'   isProband = c(1, 0)
#' )
#' 
#' # Initialize ages with random values between 20 and 94
#' result <- imputeAgesInit(data, threshold = 20, max_age = 94)
#' 
#' # Access the results
#' imputed_data <- result$data
#' missing_indices <- result$na_indices
imputeAgesInit <- function(data, threshold, max_age) {
  na_indices <- which(is.na(data$age))
  data$age[na_indices] <- runif(length(na_indices), threshold, max_age)
  return(list(data = data, na_indices = na_indices))
}

#' Calculate Empirical Age Density
#'
#' @description
#' Calculates empirical age density distributions for different subgroups in the data,
#' separated by sex and genetic testing status.
#'
#' @param data A data frame containing the family data
#' @param aff_column Name of the affection status column
#' @param age_column Name of the age column
#' @param sex_column Name of the sex column
#' @param geno_column Name of the genotype column
#' @param n_points Number of points to use in density estimation
#' @param sex_specific Logical; whether to calculate sex-specific densities
#'
#' @return A list of density objects for different subgroups (tested/untested, male/female)
#' @export
calculateEmpiricalDensity <- function(data, aff_column = "aff", age_column = "age", sex_column = "sex", 
                                    geno_column = "geno", n_points = 10000, sex_specific = TRUE) {
  
  # Helper function to check if geno is missing or empty
  is_geno_missing <- function(x) is.na(x) | x == ""
  
  # Helper function to safely calculate density
  density <- function(x, n) {
    if (length(x) < 2) return(NULL)
    tryCatch({
      density(na.omit(x), n = n)
    }, error = function(e) {
      return(NULL)
    })
  }
  
  if (sex_specific) {
    # Original sex-specific code
    empirical_density <- list(
      male_tested = density(
        subset(data, data[[aff_column]] == 0 & data[[sex_column]] == 1 & !is_geno_missing(data[[geno_column]]))[[age_column]], 
        n_points
      ),
      male_untested = density(
        subset(data, data[[aff_column]] == 0 & data[[sex_column]] == 1 & is_geno_missing(data[[geno_column]]))[[age_column]], 
        n_points
      ),
      female_tested = density(
        subset(data, data[[aff_column]] == 0 & data[[sex_column]] == 2 & !is_geno_missing(data[[geno_column]]))[[age_column]], 
        n_points
      ),
      female_untested = density(
        subset(data, data[[aff_column]] == 0 & data[[sex_column]] == 2 & is_geno_missing(data[[geno_column]]))[[age_column]], 
        n_points
      )
    )
  } else {
    # Non-sex-specific: calculate overall densities
    empirical_density <- list(
      overall_tested = density(
        subset(data, data[[aff_column]] == 0 & !is_geno_missing(data[[geno_column]]))[[age_column]], 
        n_points
      ),
      overall_untested = density(
        subset(data, data[[aff_column]] == 0 & is_geno_missing(data[[geno_column]]))[[age_column]], 
        n_points
      )
    )
  }
  
  # Handle NULL cases appropriately based on sex_specific setting
  if (sex_specific) {
    # Original NULL handling for sex-specific case
    if (is.null(empirical_density$male_tested)) empirical_density$male_tested <- empirical_density$male_untested
    if (is.null(empirical_density$male_untested)) empirical_density$male_untested <- empirical_density$male_tested
    if (is.null(empirical_density$female_tested)) empirical_density$female_tested <- empirical_density$female_untested
    if (is.null(empirical_density$female_untested)) empirical_density$female_untested <- empirical_density$female_tested
    
    # Create uniform fallback if needed
    if (is.null(empirical_density$male_tested) && is.null(empirical_density$male_untested)) {
      x <- seq(1, 100, length.out = n_points)
      y <- rep(1/length(x), length(x))
      empirical_density$male_tested <- empirical_density$male_untested <- list(x = x, y = y)
    }
    if (is.null(empirical_density$female_tested) && is.null(empirical_density$female_untested)) {
      x <- seq(1, 100, length.out = n_points)
      y <- rep(1/length(x), length(x))
      empirical_density$female_tested <- empirical_density$female_untested <- list(x = x, y = y)
    }
  } else {
    # NULL handling for non-sex-specific case
    if (is.null(empirical_density$overall_tested)) empirical_density$overall_tested <- empirical_density$overall_untested
    if (is.null(empirical_density$overall_untested)) empirical_density$overall_untested <- empirical_density$overall_tested
    
    # Create uniform fallback if needed
    if (is.null(empirical_density$overall_tested) && is.null(empirical_density$overall_untested)) {
      x <- seq(1, 100, length.out = n_points)
      y <- rep(1/length(x), length(x))
      empirical_density$overall_tested <- empirical_density$overall_untested <- list(x = x, y = y)
    }
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
#' This function draws ages using the inverse CDF method from empirical density data,
#' based on sex and whether the individual was tested.
#'
#' @param empirical_density A list of density objects containing the empirical density of ages for different groups.
#' @param sex Numeric, the sex of the individual (1 for male, 2 for female).
#' @param tested Logical, indicating whether the individual was tested (has a non-NA 'geno' value).
#' @param sex_specific Logical, indicating whether the imputation should be sex-specific. Default is TRUE.
#'
#' @return A single age value drawn from the appropriate empirical density data.
#'
drawEmpirical <- function(empirical_density, sex, tested, sex_specific = TRUE) {
  u <- runif(1)
  
  if (sex_specific) {
    # If sex is NA in sex-specific mode, return constant value
    if (is.na(sex)) return(50)
    
    density_data <- if (sex == 1) {  # Male
      if(tested) empirical_density$male_tested else empirical_density$male_untested
    } else if (sex == 2) {  # Female
      if(tested) empirical_density$female_tested else empirical_density$female_untested
    } else {
      stop("Invalid sex value: must be 1, 2, or NA")
    }
  } else {
    # Non-sex-specific: use overall densities
    density_data <- if(tested) empirical_density$overall_tested else empirical_density$overall_untested
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

#' Impute Missing Ages in Family-Based Data
#'
#' @examples
#' # Create sample data
#' data <- data.frame(
#'   family = c(1, 1, 1),
#'   individual = c(1, 2, 3),
#'   father = c(NA, 1, 1),
#'   mother = c(NA, 2, 2),
#'   sex = c(1, 2, 1),
#'   aff = c(1, 0, NA),
#'   age = c(45, NA, 20),
#'   geno = c("1/2", NA, NA),
#'   isProband = c(1, 0, 0)
#' )
#' 
#' # Initialize parameters
#' na_indices <- which(is.na(data$age))
#' max_age <- 94
#' geno_freq <- c(0.999, 0.001)
#' trans <- matrix(c(1, 0, 0.5, 0.5), nrow = 2)
#' lik <- matrix(1, nrow = nrow(data), ncol = 2)
#' 
#' # Impute ages
#' imputed_data <- imputeAges(
#'   data = data,
#'   na_indices = na_indices,
#'   max_age = max_age,
#'   geno_freq = geno_freq,
#'   trans = trans,
#'   lik = lik
#' )

#' @examples
#' # Create sample data
#' data <- data.frame(
#'   sex = c(1, 1, 2, 2),
#'   aff = c(0, 0, 0, 0),
#'   age = c(45, 50, 35, 40),
#'   geno = c("1/2", NA, "1/2", NA)
#' )
#' 
#' # Calculate density with sex-specific parameters
#' density_sex <- calculateEmpiricalDensity(data, sex_specific = TRUE)
#' 
#' # Calculate density without sex-specific parameters
#' density_nosex <- calculateEmpiricalDensity(data, sex_specific = FALSE)