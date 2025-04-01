#' Execution of a Single Chain in Metropolis-Hastings for Cancer Risk Estimation
#'
#' Performs a single chain execution in the Metropolis-Hastings algorithm for Bayesian inference,
#' specifically tailored for cancer risk estimation. This function can handle both sex-specific and
#' non-sex-specific scenarios.
#'
#' @param seed Integer, the seed for the random number generator to ensure reproducibility.
#' @param n_iter Integer, the number of iterations to perform in the Metropolis-Hastings algorithm.
#' @param burn_in Integer, the number of initial iterations to discard (burn-in period).
#' @param chain_id Integer, the identifier for the chain being executed.
#' @param ncores Integer, the number of cores to use for parallel computation.
#' @param data Data frame, containing family and genetic information used in the analysis.
#' @param twins Information on monozygous twins or triplets in the pedigrees.
#' @param max_age Integer, the maximum age considered in the analysis.
#' @param baseline_data Numeric matrix or vector, containing baseline risk estimates for different ages and sexes.
#' @param prior_distributions List, containing prior distributions for the parameters being estimated.
#' @param prev Numeric, the prevalence of the risk allele in the population.
#' @param median_max Logical, indicates if the maximum median age should be used for the Weibull distribution.
#' @param BaselineNC Logical, indicates if non-carrier penetrance should be based on SEER data.
#' @param var Numeric, the variance for the proposal distribution in the Metropolis-Hastings algorithm.
#' @param age_imputation Logical, indicates if age imputation should be performed.
#' @param imp_interval Integer, the interval at which age imputation should be performed when age_imputation = TRUE.
#' @param remove_proband Logical, indicates if the proband should be removed from the analysis.
#' @param sex_specific Logical, indicates if the analysis should differentiate by sex.
#'
#' @return A list containing samples, log likelihoods, log-acceptance ratio, and rejection rate for each iteration.
#'
#' @examples
#' # Create sample data in PanelPRO format
#' data <- data.frame(
#'   ID = 1:10,
#'   PedigreeID = rep(1, 10),
#'   Sex = c(0, 1, 0, 1, 0, 1, 0, 1, 0, 1), # 0=female, 1=male
#'   MotherID = c(NA, NA, 1, 1, 3, 3, 5, 5, 7, 7),
#'   FatherID = c(NA, NA, 2, 2, 4, 4, 6, 6, 8, 8),
#'   isProband = c(1, rep(0, 9)),
#'   CurAge = c(45, 35, 55, 40, 50, 45, 60, 38, 52, 42),
#'   isAff = c(1, 0, 1, 0, 1, 0, 1, 0, 1, 0),
#'   Age = c(40, NA, 50, NA, 45, NA, 55, NA, 48, NA),
#'   Geno = c(1, NA, 1, 0, 1, 0, NA, NA, 1, NA)
#' )
#'
#' # Transform data into required format
#' data <- transformDF(data)
#'
#' # Set parameters for the chain
#' seed <- 123
#' n_iter <- 10
#' burn_in <- 0.1 # 10% burn-in
#' chain_id <- 1
#' ncores <- 1
#' max_age <- 100
#'
#' # Create baseline data (simplified example)
#' baseline_data <- matrix(
#'   c(rep(0.005, max_age), rep(0.008, max_age)), # Increased baseline risks
#'   ncol = 2,
#'   dimnames = list(NULL, c("Male", "Female"))
#' )
#'
#' # Set prior distributions with carefully chosen bounds
#' prior_distributions <- list(
#'   prior_params = list(
#'     asymptote = list(g1 = 2, g2 = 3), # Mode around 0.4
#'     threshold = list(min = 20, max = 30), # Narrower range for threshold
#'     median = list(m1 = 3, m2 = 2), # Mode around 0.6
#'     first_quartile = list(q1 = 2, q2 = 3) # Mode around 0.4
#'   )
#' )
#'
#' # Create variance vector for all 8 parameters in sex-specific case
#' # Using very small variances for initial stability
#' var <- c(
#'   0.005, 0.005, # asymptotes (smaller variance since between 0-1)
#'   1, 1, # thresholds
#'   1, 1, # medians
#'   1, 1
#' ) # first quartiles
#'
#' # Run the chain
#' results <- mhChain(
#'   seed = seed,
#'   n_iter = n_iter,
#'   burn_in = burn_in,
#'   chain_id = chain_id,
#'   ncores = ncores,
#'   data = data,
#'   twins = NULL,
#'   max_age = max_age,
#'   baseline_data = baseline_data,
#'   prior_distributions = prior_distributions,
#'   prev = 0.05, # Increased prevalence
#'   median_max = FALSE, # Changed to FALSE for simpler median constraints
#'   BaselineNC = TRUE,
#'   var = var,
#'   age_imputation = FALSE,
#'   imp_interval = 10,
#'   remove_proband = TRUE,
#'   sex_specific = TRUE
#' )
#'
#' @export
mhChain <- function(seed, n_iter, burn_in, chain_id, ncores, data, twins, max_age, baseline_data,
                    prior_distributions, prev, median_max, BaselineNC, var,
                    age_imputation, imp_interval, remove_proband, sex_specific) {
  # Set seed for the chain
  set.seed(seed)

  # Prepare initial age imputation if enabled
  if (age_imputation) {
    # Initialize ages
    threshold <- prior_distributions$prior_params$threshold$min
    age_density <- calculateEmpiricalDensity(data)
    init_result <- imputeAgesInit(data, threshold, max_age)
    data <- init_result$data

    #  Extract the NA indices
    na_indices <- init_result$na_indices
  } else {
    # If age imputation is disabled, set unknown ages to 1 so they are disregarded in likelihood calculation
    data$age[is.na(data$age)] <- 1
  }

  # Calculate indices of the probands before removing them
  proband_indices <- which(data$isProband == 1)

  # Option to remove the proband after age imputation
  if (remove_proband) {
    # Instead of removing probands, set their affection status to NA
    data$aff[proband_indices] <- NA
  }

  # Initialize variables for sex-specific or non-specific model
  if (sex_specific) {
    # Process baseline risk data for males and females
    # Check if baseline data matches max_age
    if (nrow(baseline_data) != max_age) {
      warning(paste("Baseline data length (", nrow(baseline_data), ") does not match max_age (", max_age, "). Adjusting baseline data to match max_age.", sep=""))
      
      # If baseline data is longer than max_age, truncate it
      if (nrow(baseline_data) > max_age) {
        baseline_data <- baseline_data[1:max_age, ]
      } else {
        # If baseline data is shorter, extend it by repeating the last value
        last_male <- baseline_data[nrow(baseline_data), "Male"]
        last_female <- baseline_data[nrow(baseline_data), "Female"]
        extension <- data.frame(
          Male = rep(last_male, max_age - nrow(baseline_data)),
          Female = rep(last_female, max_age - nrow(baseline_data))
        )
        baseline_data <- rbind(baseline_data, extension)
      }
    }
    
    baseline_male <- as.numeric(baseline_data[, "Male"])
    baseline_female <- as.numeric(baseline_data[, "Female"])
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

    baseline_mid_male <- which(baseline_male_cum >= midpoint_prob_male)[1]
    baseline_mid_female <- which(baseline_female_cum >= midpoint_prob_female)[1]

    # Function to initialize parameters
    draw_initial_params <- function(data, prior_distributions) {
      data_male_affected <- data[data$sex == 1 & data$aff == 1, ]
      data_female_affected <- data[data$sex == 2 & data$aff == 1, ]

      lower_bound <- prior_distributions$prior_params$threshold$min
      upper_bound <- prior_distributions$prior_params$threshold$max

      # Add better NA handling for threshold calculations
      threshold_male <- if(length(data_male_affected$age) > 0 && sum(!is.na(data_male_affected$age)) > 0) {
        quantile(data_male_affected$age, 0.1, na.rm = TRUE)
      } else {
        (lower_bound + upper_bound) / 2  # Use middle of allowed range as default
      }
      
      threshold_female <- if(length(data_female_affected$age) > 0 && sum(!is.na(data_female_affected$age)) > 0) {
        quantile(data_female_affected$age, 0.1, na.rm = TRUE)
      } else {
        (lower_bound + upper_bound) / 2  # Use middle of allowed range as default
      }

      threshold_male <- pmax(pmin(threshold_male, upper_bound, na.rm = TRUE), lower_bound, na.rm = TRUE)
      threshold_female <- pmax(pmin(threshold_female, upper_bound, na.rm = TRUE), lower_bound, na.rm = TRUE)

      # Add better NA handling for median calculations
      median_male <- if(length(data_male_affected$age) > 0 && sum(!is.na(data_male_affected$age)) > 0) {
        median(data_male_affected$age, na.rm = TRUE)
      } else {
        50  # Default value when no data available
      }
      
      median_female <- if(length(data_female_affected$age) > 0 && sum(!is.na(data_female_affected$age)) > 0) {
        median(data_female_affected$age, na.rm = TRUE)
      } else {
        50  # Default value when no data available
      }

      # Add better NA handling for first quartile calculations
      first_quartile_male <- if(length(data_male_affected$age) > 0 && sum(!is.na(data_male_affected$age)) > 0) {
        q25 <- quantile(data_male_affected$age, probs = 0.25, na.rm = TRUE)
        max(min(q25, median_male - 1), threshold_male + 1)
      } else {
        threshold_male + 0.3 * (median_male - threshold_male)  # Default position between threshold and median
      }
      
      first_quartile_female <- if(length(data_female_affected$age) > 0 && sum(!is.na(data_female_affected$age)) > 0) {
        q25 <- quantile(data_female_affected$age, probs = 0.25, na.rm = TRUE)
        max(min(q25, median_female - 1), threshold_female + 1)
      } else {
        threshold_female + 0.3 * (median_female - threshold_female)  # Default position between threshold and median
      }
      
      # Ensure parameter relationships are maintained
      first_quartile_male <- max(first_quartile_male, threshold_male + 1)
      first_quartile_female <- max(first_quartile_female, threshold_female + 1)
      
      median_male <- max(median_male, first_quartile_male + 1)
      median_female <- max(median_female, first_quartile_female + 1)

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
  } else {
    # Use the baseline data directly as a vector for non-sex-specific
    # Check if baseline data matches max_age
    if (length(baseline_data) != max_age) {
      warning(paste("Baseline data length (", length(baseline_data), ") does not match max_age (", max_age, "). Adjusting baseline data to match max_age.", sep=""))
      
      # If baseline data is longer than max_age, truncate it
      if (length(baseline_data) > max_age) {
        baseline_data <- baseline_data[1:max_age]
      } else {
        # If baseline data is shorter, extend it by repeating the last value
        last_value <- baseline_data[length(baseline_data)]
        extension <- rep(last_value, max_age - length(baseline_data))
        baseline_data <- c(baseline_data, extension)
      }
    }
    
    baseline_cum <- cumsum(baseline_data)
    baseline_df <- data.frame(
      age = 1:length(baseline_data),
      cum_prob = baseline_cum / max(baseline_cum)
    )

    midpoint_prob <- baseline_cum[length(baseline_cum)] / 2
    baseline_mid <- which(baseline_cum >= midpoint_prob)[1]

    # Function to initialize parameters
    draw_initial_params <- function(data, prior_distributions) {
      data_affected <- data[data$aff == 1, ]

      lower_bound <- prior_distributions$prior_params$threshold$min
      upper_bound <- prior_distributions$prior_params$threshold$max

      # Add better NA handling for threshold calculation
      threshold <- if(length(data_affected$age) > 0 && sum(!is.na(data_affected$age)) > 0) {
        quantile(data_affected$age, 0.1, na.rm = TRUE)
      } else {
        (lower_bound + upper_bound) / 2  # Use middle of allowed range as default
      }

      threshold <- pmax(pmin(threshold, upper_bound, na.rm = TRUE), lower_bound, na.rm = TRUE)

      # Add better NA handling for median calculation
      median_age <- if(length(data_affected$age) > 0 && sum(!is.na(data_affected$age)) > 0) {
        median(data_affected$age, na.rm = TRUE)
      } else {
        50  # Default value when no data available
      }

      # Add better NA handling for first quartile calculation
      first_quartile <- if(length(data_affected$age) > 0 && sum(!is.na(data_affected$age)) > 0) {
        q25 <- quantile(data_affected$age, probs = 0.25, na.rm = TRUE)
        min(q25, median_age - 1)
      } else {
        threshold + 0.3 * (median_age - threshold)  # Default position between threshold and median
      }
      
      # Ensure parameter relationships are maintained 
      first_quartile <- max(first_quartile, threshold + 1)
      median_age <- max(median_age, first_quartile + 1)

      asymptote <- runif(1, max(baseline_cum), 1)

      return(list(
        asymptote = asymptote,
        threshold = threshold,
        median = median_age,
        first_quartile = first_quartile
      ))
    }
  }

  # Initialize parameters using the function draw_initial_params
  initial_params <- draw_initial_params(data = data, prior_distributions = prior_distributions)
  params_current <- initial_params
  current_states <- list()

  num_pars <- if (sex_specific) 8 else 4 # Number of parameters depends on the model type
  C <- diag(var)
  sd <- 2.38^2 / num_pars
  eps <- 0.01

  # Output lists
  if (sex_specific) {
    out <- list(
      asymptote_male_samples = numeric(n_iter),
      asymptote_female_samples = numeric(n_iter),
      threshold_male_samples = numeric(n_iter),
      threshold_female_samples = numeric(n_iter),
      median_male_samples = numeric(n_iter),
      median_female_samples = numeric(n_iter),
      first_quartile_male_samples = numeric(n_iter),
      first_quartile_female_samples = numeric(n_iter),
      loglikelihood_current = numeric(n_iter),
      loglikelihood_proposal = numeric(n_iter),
      logprior_current = numeric(n_iter),
      logprior_proposal = numeric(n_iter),
      log_acceptance_ratio = numeric(n_iter),
      rejection_rate = numeric(n_iter),
      C = vector("list", n_iter)
    )
  } else {
    out <- list(
      asymptote_samples = numeric(n_iter),
      threshold_samples = numeric(n_iter),
      median_samples = numeric(n_iter),
      first_quartile_samples = numeric(n_iter),
      loglikelihood_current = numeric(n_iter),
      loglikelihood_proposal = numeric(n_iter),
      logprior_current = numeric(n_iter),
      logprior_proposal = numeric(n_iter),
      log_acceptance_ratio = numeric(n_iter),
      rejection_rate = numeric(n_iter),
      C = vector("list", n_iter)
    )
  }

  # Function to calculate the (log) prior probabilities
  calculate_log_prior <- function(params, prior_distributions, max_age) {
    prior_params <- prior_distributions$prior_params

    if (sex_specific) {
      # Add checks to prevent division by zero
      if (max_age <= params$threshold_male) {
        return(-Inf)  # Invalid parameter, reject immediately
      }
      if (max_age <= params$threshold_female) {
        return(-Inf)  # Invalid parameter, reject immediately
      }
      if (params$median_male <= params$threshold_male) {
        return(-Inf)  # Invalid parameter, reject immediately
      }
      if (params$median_female <= params$threshold_female) {
        return(-Inf)  # Invalid parameter, reject immediately
      }
      
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
    } else {
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
    }

    return(log_prior_total)
  }

  num_rejections <- 0
  message("Starting Chain ", chain_id)

  # Initialize the model
  # geno_freq represents the frequency of the risk type and its complement in the population
  # prev is the frequency of the risk allele.
  geno_freq <- c(1 - prev, prev)

  # trans is a transition matrix that defines the probabilities of allele transmission from parents to offspring
  # We are assuming that homozygous genotype is not viable
  # Here, the rows correspond to the 4 possible joint parental genotypes and the two columns correspond to the
  # two possible offspring genotypes. Each number is the conditional probability of the offspring genotype, given the parental genotypes.
  # The first column corresponds to the wildtype and the second column to the heterozygous carrier (i.e. mutated) for the offspring.
  trans <- matrix(
    c(
      1, 0, # both parents are wild type
      0.5, 0.5, # mother is wildtype and father is a heterozygous carrier
      0.5, 0.5, # father is wildtype and mother is a heterozygous carrier
      1 / 3, 2 / 3 # both parents are heterozygous carriers
    ),
    nrow = 4, ncol = 2, byrow = TRUE
  )

  # Main loop of Metropolis-Hastings algorithm
  for (i in 1:n_iter) {
    if (sex_specific) {
      # Calculate Weibull parameters for male and female
      weibull_params_male <- calculate_weibull_parameters(params_current$median_male, params_current$first_quartile_male, params_current$threshold_male)
      weibull_params_female <- calculate_weibull_parameters(params_current$median_female, params_current$first_quartile_female, params_current$threshold_female)

      # Impute ages every imp_interval iterations, if age_imputation is TRUE.
      # Skip the first iteration since loglikelihood_current isn't available yet
      if (age_imputation && i > 1 && i %% imp_interval == 0) {
        data <- imputeAges(
          data = data,
          na_indices = na_indices,
          baseline_male = baseline_male_df,
          baseline_female = baseline_female_df,
          alpha_male = weibull_params_male$alpha,
          beta_male = weibull_params_male$beta,
          delta_male = params_current$threshold_male,
          alpha_female = weibull_params_female$alpha,
          beta_female = weibull_params_female$beta,
          delta_female = params_current$threshold_female,
          max_age = max_age,
          sex_specific = TRUE,
          geno_freq = geno_freq,
          trans = trans,
          lik = loglikelihood_current$penet
        )
      }
      # Current parameter vector for sex-specific model
      params_vector <- c(
        params_current$asymptote_male, params_current$asymptote_female,
        params_current$threshold_male, params_current$threshold_female,
        params_current$median_male, params_current$median_female,
        params_current$first_quartile_male, params_current$first_quartile_female
      )

      # Draw Proposals
      proposal_vector <- mvrnorm(1, mu = params_vector, Sigma = C)

      # Ensure the proposals for the asymptote fall within the 0 to 1 range
      proposal_vector[1] <- ifelse(proposal_vector[1] < 0, -proposal_vector[1],
        ifelse(proposal_vector[1] > 1, 2 - proposal_vector[1], proposal_vector[1])
      )
      proposal_vector[2] <- ifelse(proposal_vector[2] < 0, -proposal_vector[2],
        ifelse(proposal_vector[2] > 1, 2 - proposal_vector[2], proposal_vector[2])
      )

      # Record proposals
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

      loglikelihood_current <- mhLogLikelihood_clipp(
        params_current, data, twins, max_age,
        baseline_data, prev, geno_freq, trans, BaselineNC, ncores
      )
      logprior_current <- calculate_log_prior(params_current, prior_distributions, max_age)
    } else {
      # Non-sex-specific
      weibull_params <- calculate_weibull_parameters(params_current$median, params_current$first_quartile, params_current$threshold)

      # Impute ages only after warmup_iterations, if age_imputation is TRUE
      if (age_imputation && i %% imp_interval == 0) {
        data <- imputeAges(
          data = data,
          na_indices = na_indices,
          baseline = baseline_df,
          alpha = weibull_params$alpha,
          beta = weibull_params$beta,
          delta = params_current$threshold,
          max_age = max_age,
          sex_specific = FALSE,
          geno_freq = geno_freq,
          trans = trans,
          lik = loglikelihood_current$penet
        )
      }

      # Current parameter vector for non-sex-specific model
      params_vector <- c(params_current$asymptote, params_current$threshold, params_current$median, params_current$first_quartile)

      # Draw Proposals. Here C is 4x4
      proposal_vector <- mvrnorm(1, mu = params_vector, Sigma = C)

      # Ensure the proposals for the asymptote fall within the 0 to 1 range
      proposal_vector[1] <- ifelse(proposal_vector[1] < 0, -proposal_vector[1],
        ifelse(proposal_vector[1] > 1, 2 - proposal_vector[1], proposal_vector[1])
      )

      # Record proposals
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

      loglikelihood_current <- mhLogLikelihood_clipp_noSex(
        params_current, data, twins, max_age, baseline_data, prev, geno_freq, trans, BaselineNC, ncores
      )
      logprior_current <- calculate_log_prior(params_current, prior_distributions, max_age)
    }

    # Record the outputs of the evaluation for the current set of parameters
    out$loglikelihood_current[i] <- loglikelihood_current$loglik
    out$logprior_current[i] <- logprior_current

    # Initialize valid proposal
    valid_proposal <- TRUE

    # Explicit checks for sex-specific parameters
    if (sex_specific) {
      # Asymptote checks (male and female must be strictly between 0 and 1)
      if (proposal_vector[1] <= 0 || proposal_vector[1] >= 1) valid_proposal <- FALSE
      if (proposal_vector[2] <= 0 || proposal_vector[2] >= 1) valid_proposal <- FALSE

      # Threshold checks (male and female must be within prior bounds, strictly)
      if (proposal_vector[3] <= prior_distributions$prior_params$threshold$min ||
        proposal_vector[3] >= prior_distributions$prior_params$threshold$max) {
        valid_proposal <- FALSE
      }
      if (proposal_vector[4] <= prior_distributions$prior_params$threshold$min ||
        proposal_vector[4] >= prior_distributions$prior_params$threshold$max) {
        valid_proposal <- FALSE
      }

      # First quartile must be strictly greater than the threshold (male and female)
      if (proposal_vector[7] <= proposal_vector[3]) valid_proposal <- FALSE # First quartile male <= threshold male
      if (proposal_vector[8] <= proposal_vector[4]) valid_proposal <- FALSE # First quartile female <= threshold female

      # Median must be strictly greater than the first quartile (male and female)
      if (proposal_vector[5] <= proposal_vector[7]) valid_proposal <- FALSE # Median male <= first quartile male
      if (proposal_vector[6] <= proposal_vector[8]) valid_proposal <- FALSE # Median female <= first quartile female

      # Median should not exceed baseline midpoint or max age (for both male and female)
      if (median_max) {
        if (proposal_vector[5] >= baseline_mid_male) valid_proposal <- FALSE # Median male >= baseline midpoint
        if (proposal_vector[6] >= baseline_mid_female) valid_proposal <- FALSE # Median female >= baseline midpoint
      } else {
        if (proposal_vector[5] >= max_age) valid_proposal <- FALSE # Median male >= max age
        if (proposal_vector[6] >= max_age) valid_proposal <- FALSE # Median female >= max age
      }
    } else {
      # Non-sex-specific proposal checks
      if (proposal_vector[1] <= 0 || proposal_vector[1] >= 1) valid_proposal <- FALSE # Asymptote must be strictly between 0 and 1

      # Threshold check
      if (proposal_vector[2] <= prior_distributions$prior_params$threshold$min ||
        proposal_vector[2] >= prior_distributions$prior_params$threshold$max) {
        valid_proposal <- FALSE
      }

      # First quartile must be strictly greater than the threshold
      if (proposal_vector[4] <= proposal_vector[2]) valid_proposal <- FALSE # First quartile <= threshold

      # Median must be strictly greater than the first quartile
      if (proposal_vector[3] <= proposal_vector[4]) valid_proposal <- FALSE # Median <= first quartile

      # Median baseline check
      if (median_max) {
        if (proposal_vector[3] >= baseline_mid) valid_proposal <- FALSE # Median >= baseline midpoint
      } else {
        if (proposal_vector[3] >= max_age) valid_proposal <- FALSE # Median >= max age
      }
    }

    # If valid proposal, calculate the acceptance ratio and store
    if (valid_proposal) {
      if (sex_specific) {
        loglikelihood_proposal <- mhLogLikelihood_clipp(
          params_proposal, data, twins, max_age,
          baseline_data, prev, geno_freq, trans, BaselineNC, ncores
        )
      } else {
        loglikelihood_proposal <- mhLogLikelihood_clipp_noSex(
          params_proposal, data, twins, max_age, baseline_data, prev, geno_freq, trans, BaselineNC, ncores
        )
      }
      logprior_proposal <- calculate_log_prior(params_proposal, prior_distributions, max_age)
      log_acceptance_ratio <- (loglikelihood_proposal$loglik + logprior_proposal) - (loglikelihood_current$loglik + logprior_current)

      # Metropolis-Hastings acceptance step
      if (log(runif(1)) < log_acceptance_ratio) {
        params_current <- params_proposal
      } else {
        num_rejections <- num_rejections + 1
      }
      # Record
      out$loglikelihood_proposal[i] <- loglikelihood_proposal$loglik
      out$logprior_proposal[i] <- logprior_proposal
      out$log_acceptance_ratio[i] <- log_acceptance_ratio
    } else {
      # Proposal rejected without calculating log-likelihood
      num_rejections <- num_rejections + 1
    }

    current_states[[i]] <- params_vector

    if (i > max(burn_in * n_iter, 3)) {
      C <- sd * cov(do.call(rbind, current_states)) + eps * sd * diag(num_pars)
    }

    # Store current parameters in the output (same as before)
    if (sex_specific) {
      out$asymptote_male_samples[i] <- params_current$asymptote_male
      out$asymptote_female_samples[i] <- params_current$asymptote_female
      out$threshold_male_samples[i] <- params_current$threshold_male
      out$threshold_female_samples[i] <- params_current$threshold_female
      out$median_male_samples[i] <- params_current$median_male
      out$median_female_samples[i] <- params_current$median_female
      out$first_quartile_male_samples[i] <- params_current$first_quartile_male
      out$first_quartile_female_samples[i] <- params_current$first_quartile_female
    } else {
      out$asymptote_samples[i] <- params_current$asymptote
      out$threshold_samples[i] <- params_current$threshold
      out$median_samples[i] <- params_current$median
      out$first_quartile_samples[i] <- params_current$first_quartile
    }

    out$C[[i]] <- C
  }

  out$rejection_rate <- num_rejections / n_iter

  # Return both the main results
  return(out)
}