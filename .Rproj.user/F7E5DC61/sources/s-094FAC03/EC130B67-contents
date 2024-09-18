#' Execution of a Single Chain in Metropolis-Hastings for Cancer Risk Estimation with Sex Differentiation
#'
#' Performs a single chain execution in the Metropolis-Hastings algorithm for Bayesian inference,
#' specifically tailored for cancer risk estimation. It estimates parameters related to cancer penetrance
#' based on family data, genetic information, and baseline database estimates.
#'
#' @param seed Integer, the seed for the random number generator to ensure reproducibility.
#' @param n_iter Integer, the number of iterations to perform in the Metropolis-Hastings algorithm.
#' @param burn_in Integer, the number of initial iterations to discard (burn-in period).
#' @param chain_id Integer, the identifier for the chain being executed.
#' @param ncores Integer, the number of cores to use for parallel computation.
#' @param data Data frame, containing family and genetic information used in the analysis.
#' @param twins Information on monozygous twins or triplets in the pedigrees.
#' @param max_age Integer, the maximum age considered in the analysis.
#' @param baseline_data Numeric matrix, containing baseline risk estimates for different ages and sexes.
#' @param prior_distributions List, containing prior distributions for the parameters being estimated.
#' @param af Numeric, the allele frequency of the risk allele in the population.
#' @param median_max Numeric, the maximum median age for the Weibull distribution.
#' @param max_penetrance Numeric, the maximum penetrance value allowed.
#' @param BaselineNC Logical, indicates if non-carrier penetrance should be based on SEER data.
#' @param var Numeric, the variance for the proposal distribution in the Metropolis-Hastings algorithm.
#' @param ageImputation Logical, indicates if age imputation should be performed.
#' @param removeProband Logical, indicates if the proband should be removed from the analysis.
#'
#' @return A list containing samples, log likelihoods, acceptance ratio, and rejection rate for each iteration.
#' 
mhChain <- function(seed, n_iter, burn_in, chain_id, ncores, data, twins, max_age, baseline_data,
                    prior_distributions, af, median_max, max_penetrance, BaselineNC, var,
                    ageImputation, removeProband) {
  
  # Set seed for the chain
  set.seed(seed)
  
  # Calculate empirical age density for affected individuals
  age_density <- calculateEmpiricalDensity(data, aff_column = "aff", age_column = "age")
  
  # Prepare initial age imputation if enabled
  if (ageImputation) {
    data <- calcPedDegree(data)
    # Initialize ages
    threshold <- prior_distributions$prior_params$threshold$min
    init_result <- imputeAgesInit(data, threshold, max_age)
    data <- init_result$data
    na_indices <- init_result$na_indices
  } else {
    # If age imputation is disabled, set unknown ages to 1 so they are disregarded in likelihood calculation
    data$age[is.na(data$age)] <- 1
  }
  
  # Calculate indices of the probands before removing them
  proband_indices <- which(data$isProband == 1)
  
  # Option to remove the proband after age imputation
  if (removeProband) {
    
    # Sort proband_indices in descending order
    # This ensures that removing rows doesn't affect the remaining indices
    proband_indices <- sort(proband_indices, decreasing = TRUE)
    
    # Remove rows corresponding to the proband indices
    data <- data[-proband_indices, ]
    
    # Adjust na_indices
    if (exists("na_indices")) {
      original_na_indices <- na_indices
      
      for (proband_index in proband_indices) {
        # Remove the current proband index from na_indices
        original_na_indices <- original_na_indices[original_na_indices != proband_index]
        
        # Decrement all indices greater than the current proband_index
        original_na_indices[original_na_indices > proband_index] <- original_na_indices[original_na_indices > proband_index] - 1
      }
      
      # Update na_indices with the adjusted indices
      na_indices <- original_na_indices
    }
  }
  
  # Process baseline risk data
  baseline_male <- as.numeric(baseline_data[,"Male"])
  baseline_female <- as.numeric(baseline_data[,"Female"])
  
  baseline_male_cum <- cumsum(baseline_male)
  baseline_female_cum <- cumsum(baseline_female)
  
  baseline_male_df <- data.frame(
    age = 1:length(baseline_male),
    cum_prob = baseline_male_cum / max(baseline_male_cum)
  )
  baseline_female_df <- data.frame(
    age = 1:length(baseline_female),
    cum_prob = baseline_female_cum / max(baseline_female_cum)
  )
  
  midpoint_prob_male <- baseline_male_cum[length(baseline_male_cum)] / 2
  midpoint_prob_female <- baseline_female_cum[length(baseline_female_cum)] / 2
  
  midpoint_index_male <- which(baseline_male_cum >= midpoint_prob_male)[1]
  midpoint_index_female <- which(baseline_female_cum >= midpoint_prob_female)[1]
  
  baseline_mid_male <- midpoint_index_male
  baseline_mid_female <- midpoint_index_female
  
  # Function to initialize the Weibull parameters using empirical data
  draw_initial_params <- function(data, prior_distributions) {

    # Filter data by sex and affected status
    data_male_affected <- data[data$sex == 1 & data$aff == 1, ]
    data_female_affected <- data[data$sex == 2 & data$aff == 1, ]
    
    # Calculate threshold (minimum age), median, and first quartile by sex among affected individuals
    lower_bound <- prior_distributions$prior_params$threshold$min
    upper_bound <- prior_distributions$prior_params$threshold$max
    
    # Initialize using the first decile for the threshold 
    threshold_male <- ifelse(length(data_male_affected$age) > 0,
                             quantile(data_male_affected$age, 0.1, na.rm = TRUE), NA)
    
    threshold_female <- ifelse(length(data_female_affected$age) > 0,
                               quantile(data_female_affected$age, 0.1, na.rm = TRUE), NA)
    
    threshold_male <- pmax(pmin(threshold_male, upper_bound, na.rm = TRUE), lower_bound, na.rm = TRUE)
    threshold_female <- pmax(pmin(threshold_female, upper_bound, na.rm = TRUE), lower_bound, na.rm = TRUE)
    
    median_male <- ifelse(length(data_male_affected$age) > 0,
                          median(data_male_affected$age, na.rm = TRUE), 50
    )
    median_female <- ifelse(length(data_female_affected$age) > 0,
                            median(data_female_affected$age, na.rm = TRUE), 50
    )
    
    first_quartile_male <- ifelse(length(data_male_affected$age) > 0,
                                  min(quantile(data_male_affected$age, probs = 0.25, na.rm = TRUE), median_male - 1), 40
    )
    first_quartile_female <- ifelse(length(data_female_affected$age) > 0,
                                    min(quantile(data_female_affected$age, probs = 0.25, na.rm = TRUE), median_female - 1), 40
    )
    
    asymptote_male <- runif(1, max(baseline_male_cum), 1)
    asymptote_female <- runif(1, max(baseline_female_cum), 1)
    
    return(list(
      asymptote_male = asymptote_male,
      asymptote_female = asymptote_female,
      threshold_male = threshold_male,
      threshold_female = threshold_female,
      median_male = median_male,
      median_female = median_female,
      first_quartile_male = first_quartile_male,
      first_quartile_female = first_quartile_female
    ))
  }
  
  # Initialize Parameters using the function draw_initial_params
  initial_params <- draw_initial_params(data = data, prior_distributions = prior_distributions)
  params_current <- initial_params
  current_states <- list()
  
  num_pars <- 8
  C <- diag(var)
  sd <- 2.38^2 / num_pars
  eps <- 0.01
  
  out <- list(
    asymptote_male_samples = numeric(n_iter),
    asymptote_female_samples = numeric(n_iter),
    asymptote_male_proposals = numeric(n_iter),
    asymptote_female_proposals = numeric(n_iter),
    threshold_male_samples = numeric(n_iter),
    threshold_female_samples = numeric(n_iter),
    threshold_male_proposals = numeric(n_iter),
    threshold_female_proposals = numeric(n_iter),
    median_male_samples = numeric(n_iter),
    median_female_samples = numeric(n_iter),
    median_male_proposals = numeric(n_iter),
    median_female_proposals = numeric(n_iter),
    first_quartile_male_samples = numeric(n_iter),
    first_quartile_female_samples = numeric(n_iter),
    first_quartile_male_proposals = numeric(n_iter),
    first_quartile_female_proposals = numeric(n_iter),
    loglikelihood_current = numeric(n_iter),
    loglikelihood_proposal = numeric(n_iter),
    logprior_current = numeric(n_iter),
    logprior_proposal = numeric(n_iter),
    acceptance_ratio = numeric(n_iter),
    rejection_rate = numeric(n_iter),
    C = vector("list", n_iter),
    data = vector("list", n_iter)
  )
  
  num_rejections <- 0
  cat("Starting Chain", chain_id, "\n")
  
  # Function to calculate the (log) prior probabilities using the fixed prior distributions and specified prior parameters. 
  calculate_log_prior <- function(params, prior_distributions, max_age) {
    prior_params <- prior_distributions$prior_params
    
    scaled_asymptote_male <- params$asymptote_male
    scaled_asymptote_female <- params$asymptote_female
    
    scaled_threshold_male <- params$threshold_male
    scaled_threshold_female <- params$threshold_female
    
    scaled_median_male <- (params$median_male - params$threshold_male) / (max_age - params$threshold_male)
    scaled_median_female <- (params$median_female - params$threshold_female) / (max_age - params$threshold_female)
    
    scaled_first_quartile_male <- (params$first_quartile_male - params$threshold_male) /
      (params$median_male - params$threshold_male)
    scaled_first_quartile_female <- (params$first_quartile_female - params$threshold_female) /
      (params$median_female - params$threshold_female)
    
    log_prior_asymptote_male <- dbeta(scaled_asymptote_male, prior_params$asymptote$g1, prior_params$asymptote$g2, log = TRUE)
    log_prior_asymptote_female <- dbeta(scaled_asymptote_female, prior_params$asymptote$g1, prior_params$asymptote$g2, log = TRUE)
    
    log_prior_threshold_male <- dunif(scaled_threshold_male, prior_params$threshold$min, prior_params$threshold$max, log = TRUE)
    log_prior_threshold_female <- dunif(scaled_threshold_female, prior_params$threshold$min, prior_params$threshold$max, log = TRUE)
    
    log_prior_median_male <- dbeta(scaled_median_male, prior_params$median$m1, prior_params$median$m2, log = TRUE)
    log_prior_median_female <- dbeta(scaled_median_female, prior_params$median$m1, prior_params$median$m2, log = TRUE)
    
    log_prior_first_quartile_male <- dbeta(scaled_first_quartile_male, prior_params$first_quartile$q1, prior_params$first_quartile$q2, log = TRUE)
    log_prior_first_quartile_female <- dbeta(scaled_first_quartile_female, prior_params$first_quartile$q1, prior_params$first_quartile$q2, log = TRUE)
    
    log_prior_total <- log_prior_asymptote_male + log_prior_asymptote_female +
      log_prior_threshold_male + log_prior_threshold_female +
      log_prior_median_male + log_prior_median_female +
      log_prior_first_quartile_male + log_prior_first_quartile_female
    
    return(log_prior_total)
  }
  
  # Run n_iter iterations of the adaptive Metropolis-Hastings algorithm
  for (i in 1:n_iter) {

      # Calculate Weibull parameters from current parameters
      weibull_params_male <- calculate_weibull_parameters(params_current$median_male, params_current$first_quartile_male, params_current$threshold_male)
      alpha_male <- weibull_params_male$alpha
      beta_male <- weibull_params_male$beta
      delta_male <- params_current$threshold_male
      
      weibull_params_female <- calculate_weibull_parameters(params_current$median_female, params_current$first_quartile_female, params_current$threshold_female)
      alpha_female <- weibull_params_female$alpha
      beta_female <- weibull_params_female$beta
      delta_female <- params_current$threshold_female
      
      # Impute ages at each iteration based on current parameters
      if (ageImputation) {
        data <- imputeAges(
          data = data, 
          na_indices = na_indices, 
          baseline_male = baseline_male_df, 
          baseline_female = baseline_female_df,
          alpha_male = alpha_male, 
          beta_male = beta_male, 
          delta_male = delta_male,
          alpha_female = alpha_female, 
          beta_female = beta_female, 
          delta_female = delta_female,
          empirical_density = age_density, 
          max_age = max_age, 
          sex_specific = TRUE
        )
      }

      # Store the current parameter values
      params_vector <- c(
        params_current$asymptote_male, params_current$asymptote_female,
        params_current$threshold_male, params_current$threshold_female,
        params_current$median_male, params_current$median_female,
        params_current$first_quartile_male, params_current$first_quartile_female
      )
      
      # Generate the proposal parameters from a multivariate normal distribution centered around current parameters
      proposal_vector <- mvrnorm(1, mu = params_vector, Sigma = C)
      
      # Ensure the proposals for the asymptote fall within the 0 to 1 range
      proposal_vector[1] <- ifelse(proposal_vector[1] < 0, -proposal_vector[1],
                                   ifelse(proposal_vector[1] > 1, 2 - proposal_vector[1], proposal_vector[1])
      )
      proposal_vector[2] <- ifelse(proposal_vector[2] < 0, -proposal_vector[2],
                                   ifelse(proposal_vector[2] > 1, 2 - proposal_vector[2], proposal_vector[2])
      )
      
      out$asymptote_male_proposals[i] <- proposal_vector[1]
      out$asymptote_female_proposals[i] <- proposal_vector[2]
      out$threshold_male_proposals[i] <- proposal_vector[3]
      out$threshold_female_proposals[i] <- proposal_vector[4]
      out$median_male_proposals[i] <- proposal_vector[5]
      out$median_female_proposals[i] <- proposal_vector[6]
      out$first_quartile_male_proposals[i] <- proposal_vector[7]
      out$first_quartile_female_proposals[i] <- proposal_vector[8]
      
      params_proposal <- list(
        asymptote_male = proposal_vector[1],
        asymptote_female = proposal_vector[2],
        threshold_male = proposal_vector[3],
        threshold_female = proposal_vector[4],
        median_male = proposal_vector[5],
        median_female = proposal_vector[6],
        first_quartile_male = proposal_vector[7],
        first_quartile_female = proposal_vector[8]
      )
      
      # Evaluate the current set of parameters
      loglikelihood_current <- mhLogLikelihood_clipp(
        params_current, data, twins, max_age,
        baseline_data, af, BaselineNC, ncores
      )
      logprior_current <- calculate_log_prior(params_current, prior_distributions, max_age)
      # Record the outputs of the evaluation for the current set of parameters
      out$loglikelihood_current[i] <- loglikelihood_current
      out$logprior_current[i] <- logprior_current
      
      out$loglikelihood_proposal[i] <- NA
      out$logprior_proposal[i] <- NA
      out$acceptance_ratio[i] <- NA
      
      # Check that the proposed parameters satisfy the basic requirements
      is_rejected <- FALSE
      if (
        is.na(proposal_vector[1]) || proposal_vector[1] < 0 || proposal_vector[1] > 1 ||
        is.na(proposal_vector[2]) || proposal_vector[2] < 0 || proposal_vector[2] > 1 ||
        is.na(proposal_vector[3]) || proposal_vector[3] < 0 || proposal_vector[3] > 100 || 
        proposal_vector[3] < prior_distributions$prior_params$threshold$min || proposal_vector[3] > prior_distributions$prior_params$threshold$max ||
        is.na(proposal_vector[4]) || proposal_vector[4] < 0 || proposal_vector[4] > 100 ||
        proposal_vector[4] < prior_distributions$prior_params$threshold$min || proposal_vector[4] > prior_distributions$prior_params$threshold$max ||
        is.na(proposal_vector[5]) || proposal_vector[5] < proposal_vector[7] ||
        is.na(proposal_vector[6]) || proposal_vector[6] < proposal_vector[8] ||
        is.na(proposal_vector[7]) || proposal_vector[7] < proposal_vector[3] || proposal_vector[7] > proposal_vector[5] ||
        is.na(proposal_vector[8]) || proposal_vector[8] < proposal_vector[4] || proposal_vector[8] > proposal_vector[6] ||
        (median_max && proposal_vector[5] > baseline_mid_male) ||
        (!median_max && proposal_vector[5] > max_age) ||
        (median_max && proposal_vector[6] > baseline_mid_female) ||
        (!median_max && proposal_vector[6] > max_age)
      ) {
        is_rejected <- TRUE
        num_rejections <- num_rejections + 1
      } else {
        loglikelihood_proposal <- mhLogLikelihood_clipp(
          params_proposal, data, twins, max_age,
          baseline_data, af, BaselineNC, ncores
        )
        logprior_proposal <- calculate_log_prior(params_proposal, prior_distributions, max_age)
        
        log_acceptance_ratio <- (loglikelihood_proposal + logprior_proposal) - (loglikelihood_current + logprior_current)
        
        if (log(runif(1)) < log_acceptance_ratio) {
          params_current <- params_proposal
        } else {
          num_rejections <- num_rejections + 1
        }
        
        out$loglikelihood_proposal[i] <- loglikelihood_proposal
        out$logprior_proposal[i] <- logprior_proposal
        out$acceptance_ratio[i] <- log_acceptance_ratio
      }
      
      current_states[[i]] <- c(
        params_current$asymptote_male, params_current$asymptote_female,
        params_current$threshold_male, params_current$threshold_female,
        params_current$median_male, params_current$median_female,
        params_current$first_quartile_male, params_current$first_quartile_female
      )
      
      if (i > max(burn_in * n_iter, 3)) {
        C <- sd * cov(do.call(rbind, current_states)) + eps * sd * diag(num_pars)
      }
      
      out$asymptote_male_samples[i] <- params_current$asymptote_male
      out$asymptote_female_samples[i] <- params_current$asymptote_female
      out$threshold_male_samples[i] <- params_current$threshold_male
      out$threshold_female_samples[i] <- params_current$threshold_female
      out$median_male_samples[i] <- params_current$median_male
      out$median_female_samples[i] <- params_current$median_female
      out$first_quartile_male_samples[i] <- params_current$first_quartile_male
      out$first_quartile_female_samples[i] <- params_current$first_quartile_female
      out$C[[i]] <- C

  }
  
  out$rejection_rate <- num_rejections / n_iter
  return(out)
}

#' Execution of a Single Chain in Metropolis-Hastings for Cancer Risk Estimation without Sex Differentiation
#'
#' Performs a single chain execution in the Metropolis-Hastings algorithm for Bayesian inference,
#' specifically tailored for cancer risk estimation without considering sex differentiation. 
#' It estimates parameters related to cancer penetrance based on family data, genetic information, 
#' and baseline database estimates.
#'
#' @param seed Integer, the seed for the random number generator to ensure reproducibility.
#' @param n_iter Integer, the number of iterations to perform in the Metropolis-Hastings algorithm.
#' @param burn_in Integer, the number of initial iterations to discard (burn-in period).
#' @param chain_id Integer, the identifier for the chain being executed.
#' @param ncores Integer, the number of cores to use for parallel computation.
#' @param data Data frame, containing family and genetic information used in the analysis.
#' @param twins Information on monozygous twins or triplets in the pedigrees.
#' @param max_age Integer, the maximum age considered in the analysis.
#' @param baseline_data Numeric vector, containing baseline risk estimates for different ages.
#' @param prior_distributions List, containing prior distributions for the parameters being estimated.
#' @param af Numeric, the allele frequency of the risk allele in the population.
#' @param median_max Logical, indicates if the maximum median age should be used for the Weibull distribution.
#' @param max_penetrance Numeric, the maximum penetrance value allowed.
#' @param BaselineNC Logical, indicates if non-carrier penetrance should be based on SEER data.
#' @param var Numeric, the variance for the proposal distribution in the Metropolis-Hastings algorithm.
#' @param ageImputation Logical, indicates if age imputation should be performed.
#' @param removeProband Logical, indicates if the proband should be removed from the analysis.
#'
#' @return A list containing samples, log likelihoods, acceptance ratio, and rejection rate for each iteration.
#'
#' @export
#' 
mhChain_noSex <- function(seed, n_iter, burn_in, chain_id, ncores, data, twins, max_age, baseline_data,
                          prior_distributions, af, median_max, max_penetrance, BaselineNC, var,
                          ageImputation, removeProband) {
  # Set seed for the chain
  set.seed(seed)
  
  # Calculate empirical age density for affected individuals
  age_density <- calculateEmpiricalDensity(data, aff_column = "aff", age_column = "age")
  
  # Prepare initial age imputation if enabled
  if (ageImputation) {
    data <- calcPedDegree(data)
    # Initialize ages
    threshold <- prior_distributions$prior_params$threshold$min
    init_result <- imputeAgesInit(data, threshold, max_age)
    data <- init_result$data
    na_indices <- init_result$na_indices
  } else {
    # If age imputation is disabled, set unknown ages to 1 so they are disregarded in likelihood calculation
    data$age[is.na(data$age)] <- 1
  }
  
  # Calculate indices of the probands before removing them
  proband_indices <- which(data$isProband == 1)
  
  # Option to remove the proband after age imputation
  if (removeProband) {
    proband_indices <- sort(proband_indices, decreasing = TRUE)
    data <- data[-proband_indices, ]
    
    # Adjust na_indices
    if (exists("na_indices")) {
      original_na_indices <- na_indices
      for (proband_index in proband_indices) {
        original_na_indices <- original_na_indices[original_na_indices != proband_index]
        original_na_indices[original_na_indices > proband_index] <- original_na_indices[original_na_indices > proband_index] - 1
      }
      na_indices <- original_na_indices
    }
  }
  
  # Use the baseline data directly as a vector
  baseline_cum <- cumsum(baseline_data)
  baseline_df <- data.frame(
    age = 1:length(baseline_data),
    cum_prob = baseline_cum / max(baseline_cum)
  )
  
  midpoint_prob <- baseline_cum[length(baseline_cum)] / 2
  midpoint_index <- which(baseline_cum >= midpoint_prob)[1]
  baseline_mid <- midpoint_index
  
  # Function to initialize the Weibull parameters using empirical data
  draw_initial_params_single <- function(data, prior_distributions) {
    # Filter data by affected status
    data_affected <- data[data$aff == 1, ]
    
    # Calculate threshold (minimum age), median, and first quartile among affected individuals
    lower_bound <- prior_distributions$prior_params$threshold$min
    upper_bound <- prior_distributions$prior_params$threshold$max
    
    # Initialize using the first decile for the threshold
    threshold <- ifelse(length(data_affected$age) > 0,
                        quantile(data_affected$age, 0.1, na.rm = TRUE), NA)
    
    threshold <- pmax(pmin(threshold, upper_bound, na.rm = TRUE), lower_bound, na.rm = TRUE)
    
    median_age <- ifelse(length(data_affected$age) > 0,
                         median(data_affected$age, na.rm = TRUE), 50)
    
    first_quartile <- ifelse(length(data_affected$age) > 0,
                             min(quantile(data_affected$age, probs = 0.25, na.rm = TRUE), median_age - 1), 40)
    
    asymptote <- runif(1, max(baseline_cum), 1)
    
    return(list(
      asymptote = asymptote,
      threshold = threshold,
      median = median_age,
      first_quartile = first_quartile
    ))
  }
  
  # Initialize Parameters using the function draw_initial_params_single
  initial_params <- draw_initial_params_single(data = data, prior_distributions = prior_distributions)
  params_current <- initial_params
  current_states <- list()
  
  num_pars <- 4  # Only 4 parameters for the non-sex-specific model
  C <- diag(var)
  sd <- 2.38^2 / num_pars
  eps <- 0.01
  
  out <- list(
    asymptote_samples = numeric(n_iter),
    threshold_samples = numeric(n_iter),
    median_samples = numeric(n_iter),
    first_quartile_samples = numeric(n_iter),
    asymptote_proposals = numeric(n_iter),
    threshold_proposals = numeric(n_iter),
    median_proposals = numeric(n_iter),
    first_quartile_proposals = numeric(n_iter),
    loglikelihood_current = numeric(n_iter),
    loglikelihood_proposal = numeric(n_iter),
    logprior_current = numeric(n_iter),
    logprior_proposal = numeric(n_iter),
    acceptance_ratio = numeric(n_iter),
    rejection_rate = numeric(n_iter),
    C = vector("list", n_iter),
    data = vector("list", n_iter)
  )
  
  num_rejections <- 0
  cat("Starting Chain", chain_id, "\n")
  
  # Function to calculate the (log) prior probabilities using the fixed prior distributions and specified prior parameters.
  calculate_log_prior <- function(params, prior_distributions, max_age) {
    prior_params <- prior_distributions$prior_params
    
    scaled_asymptote <- params$asymptote
    scaled_threshold <- params$threshold
    scaled_median <- (params$median - params$threshold) / (max_age - params$threshold)
    scaled_first_quartile <- (params$first_quartile - params$threshold) /
      (params$median - params$threshold)
    
    log_prior_asymptote <- dbeta(scaled_asymptote, prior_params$asymptote$g1, prior_params$asymptote$g2, log = TRUE)
    log_prior_threshold <- dunif(scaled_threshold, prior_params$threshold$min, prior_params$threshold$max, log = TRUE)
    log_prior_median <- dbeta(scaled_median, prior_params$median$m1, prior_params$median$m2, log = TRUE)
    log_prior_first_quartile <- dbeta(scaled_first_quartile, prior_params$first_quartile$q1, prior_params$first_quartile$q2, log = TRUE)
    
    log_prior_total <- log_prior_asymptote + log_prior_threshold + log_prior_median + log_prior_first_quartile
    
    return(log_prior_total)
  }
  
  # Run n_iter iterations of the adaptive Metropolis-Hastings algorithm
  for (i in 1:n_iter) {
    # Calculate Weibull parameters from current parameters
    weibull_params <- calculate_weibull_parameters(params_current$median, params_current$first_quartile, params_current$threshold)
    alpha <- weibull_params$alpha
    beta <- weibull_params$beta
    delta <- params_current$threshold
    
    # Impute ages at each iteration based on current parameters
    if (ageImputation) {
      data <- imputeAges(
        data = data, 
        na_indices = na_indices, 
        baseline = baseline_df, 
        alpha = alpha, 
        beta = beta, 
        delta = delta,
        empirical_density = age_density, 
        max_age = max_age, 
        sex_specific = FALSE
      )
    }
    
    # Store the current parameter values
    params_vector <- c(params_current$asymptote, params_current$threshold, params_current$median, params_current$first_quartile)
    
    # Generate the proposal parameters from a multivariate normal distribution centered around current parameters
    proposal_vector <- mvrnorm(1, mu = params_vector, Sigma = C)
    
    # Ensure the proposals for the asymptote fall within the 0 to 1 range
    proposal_vector[1] <- ifelse(proposal_vector[1] < 0, -proposal_vector[1],
                                 ifelse(proposal_vector[1] > 1, 2 - proposal_vector[1], proposal_vector[1]))
    
    out$asymptote_proposals[i] <- proposal_vector[1]
    out$threshold_proposals[i] <- proposal_vector[2]
    out$median_proposals[i] <- proposal_vector[3]
    out$first_quartile_proposals[i] <- proposal_vector[4]
    
    params_proposal <- list(
      asymptote = proposal_vector[1],
      threshold = proposal_vector[2],
      median = proposal_vector[3],
      first_quartile = proposal_vector[4]
    )
    
    # Evaluate the current set of parameters
    loglikelihood_current <- mhLogLikelihood_clipp_noSex(
      params_current, data, twins, max_age, baseline_data, af, BaselineNC, ncores
    )
    logprior_current <- calculate_log_prior(params_current, prior_distributions, max_age)
    
    # Record the outputs of the evaluation for the current set of parameters
    out$loglikelihood_current[i] <- loglikelihood_current
    out$logprior_current[i] <- logprior_current
    
    out$loglikelihood_proposal[i] <- NA
    out$logprior_proposal[i] <- NA
    out$acceptance_ratio[i] <- NA
    
    # Check that the proposed parameters satisfy the basic requirements
    is_rejected <- FALSE
    if (
      is.na(proposal_vector[1]) || proposal_vector[1] < 0 || proposal_vector[1] > 1 ||
      is.na(proposal_vector[2]) || proposal_vector[2] < 0 || proposal_vector[2] > 100 ||
      proposal_vector[2] < prior_distributions$prior_params$threshold$min || proposal_vector[2] > prior_distributions$prior_params$threshold$max ||
      is.na(proposal_vector[3]) || proposal_vector[3] < proposal_vector[4] ||
      is.na(proposal_vector[4]) || proposal_vector[4] < proposal_vector[2] || proposal_vector[4] > proposal_vector[3] ||
      (median_max && proposal_vector[3] > baseline_mid) ||
      (!median_max && proposal_vector[3] > max_age)
    ) {
      is_rejected <- TRUE
      num_rejections <- num_rejections + 1
    } else {
      loglikelihood_proposal <- mhLogLikelihood_clipp_noSex(
        params_proposal, data, twins, max_age, baseline_data, af, BaselineNC, ncores
      )
      logprior_proposal <- calculate_log_prior(params_proposal, prior_distributions, max_age)
      
      log_acceptance_ratio <- (loglikelihood_proposal + logprior_proposal) - (loglikelihood_current + logprior_current)
      
      if (log(runif(1)) < log_acceptance_ratio) {
        params_current <- params_proposal
      } else {
        num_rejections <- num_rejections + 1
      }
      
      out$loglikelihood_proposal[i] <- loglikelihood_proposal
      out$logprior_proposal[i] <- logprior_proposal
      out$acceptance_ratio[i] <- log_acceptance_ratio
    }
    
    current_states[[i]] <- c(params_current$asymptote, params_current$threshold, params_current$median, params_current$first_quartile)
    
    if (i > max(burn_in * n_iter, 3)) {
      C <- sd * cov(do.call(rbind, current_states)) + eps * sd * diag(num_pars)
    }
    
    out$asymptote_samples[i] <- params_current$asymptote
    out$threshold_samples[i] <- params_current$threshold
    out$median_samples[i] <- params_current$median
    out$first_quartile_samples[i] <- params_current$first_quartile
    out$C[[i]] <- C
  }
  
  out$rejection_rate <- num_rejections / n_iter
  return(out)
}