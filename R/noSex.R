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
        data, na_indices, baseline_df, alpha, beta, delta,
        age_density, max_age
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

mhLogLikelihood_clipp_noSex <- function(paras, families, twins, max_age, baseline_data, af, BaselineNC, ncores) {
  # Extract parameters
  paras <- unlist(paras)
  gamma <- paras[1]  # Asymptote
  delta <- paras[2]  # Threshold
  given_median <- paras[3]
  given_first_quartile <- paras[4]
  
  # Calculate Weibull parameters
  params <- calculate_weibull_parameters(given_median, given_first_quartile, delta)
  alpha <- params$alpha
  beta <- params$beta
  
  # Initialize the model
  geno_freq <- c(1 - af, af)
  trans <- matrix(
    c(
      1, 0,
      0.5, 0.5,
      0.5, 0.5,
      1 / 3, 2 / 3
    ),
    nrow = 4, ncol = 2, byrow = TRUE
  )
  
  # Use the baselineRisk vector directly
  baselineRisk <- baseline_data
  
  # Calculate penetrance
  lik <- t(sapply(1:nrow(families), function(i) {
    lik_noSex(i, families, alpha, beta, delta, gamma, max_age, baselineRisk, BaselineNC)
  }))
  
  # Compute log-likelihood
  loglik <- pedigree_loglikelihood(dat = families, geno_freq = geno_freq, trans = trans, penet = lik, monozyg = twins, ncores = ncores)
  
  # Handle -Inf values
  if (is.infinite(loglik) && loglik == -Inf) {
    loglik <- -50000
  }
  
  return(loglik)
}

lik_noSex <- function(i, data, alpha, beta, delta, gamma, max_age, baselineRisk, BaselineNC) {
  if (data$age[i] == 0 || data$age[i] == 1) {
    lik.i <- c(1, 1)  # Assuming people aged 0 or 1 years are all unaffected
  } else {
    # Ensure age is within the valid range
    age_index <- min(max_age, data$age[i])
    
    # Weibull parameters for penetrance, using a single set of parameters
    survival_prob <- 1 - pweibull(max(age_index - delta, 1), shape = alpha, scale = beta) * gamma
    c.pen <- (pweibull(max(age_index - delta, 1), shape = alpha, scale = beta)
              - pweibull(max(age_index - 1 - delta, 1), shape = alpha, scale = beta)) * gamma
    
    # Extract the corresponding baseline risk for the age
    SEER_baseline_max <- baselineRisk[1:age_index]
    SEER_baseline_cum <- cumsum(baselineRisk)[age_index]
    SEER_baseline_i <- baselineRisk[age_index]
    
    # Calculate cumulative risk for non-carriers based on SEER data or other model
    if (BaselineNC == TRUE) {
      nc.pen <- SEER_baseline_i
      nc.pen.c <- prod(1 - SEER_baseline_i)
    } else {
      nc_pen_results <- calculateNCPen(
        SEER_baseline = SEER_baseline_max, alpha = alpha,
        beta = beta, delta = delta, gamma = gamma, max_age = max_age
      )
      nc.pen <- nc_pen_results$weightedCarrierRisk[age_index]
      nc.pen.c <- nc_pen_results$cumulativeProb[age_index]
    }
    
    # Penetrance calculations based on genotype and affection status
    lik.i <- c(nc.pen.c, survival_prob)  # For censored observations
    if (data$aff[i] == 1) lik.i <- c(nc.pen * nc.pen.c, c.pen)  # For affected observations
  }
  
  # Adjustment for observed genotypes
  if (data$geno[i] == "1/1") lik.i[-1] <- 0
  if (data$geno[i] == "1/2") lik.i[-2] <- 0
  
  return(lik.i)
}


####### OUTPUT FUNCTIONS #######

generate_density_plots <- function(data) {
  # Set the plotting parameters
  par(mfrow = c(3, 2), las = 1, mar = c(5, 4, 4, 2) + 0.1)
  
  # Determine which set of parameters to plot: sex-specific or non-sex-specific
  if (!is.null(data$median_male_results) || !is.null(data$median_female_results)) {
    # Plot sex-specific parameters
    plot_names <- list(
      "median_male_results" = data$median_male_results,
      "first_quartile_male_results" = data$first_quartile_male_results,
      "asymptote_male_results" = data$asymptote_male_results,
      "threshold_male_results" = data$threshold_male_results,
      "median_female_results" = data$median_female_results,
      "first_quartile_female_results" = data$first_quartile_female_results,
      "asymptote_female_results" = data$asymptote_female_results,
      "threshold_female_results" = data$threshold_female_results
    )
  } else {
    # Plot non-sex-specific parameters
    plot_names <- list(
      "median_results" = data$median_samples,
      "first_quartile_results" = data$first_quartile_samples,
      "asymptote_results" = data$asymptote_samples,
      "threshold_results" = data$threshold_samples
    )
  }
  
  # Plot each parameter that has data
  for (name in names(plot_names)) {
    param_data <- plot_names[[name]]
    
    if (is.null(param_data) || length(param_data) == 0) {
      next # Skip this iteration if the data is empty
    }
    
    mod_name <- gsub("_", " ", name)
    mod_name <- paste0(toupper(substring(mod_name, 1, 1)), substring(mod_name, 2))
    
    # Set xlim based on the name of the vector
    xlim <- if (name %in% c("median_male_results", "first_quartile_male_results", 
                            "threshold_male_results", "median_female_results", 
                            "first_quartile_female_results", 
                            "threshold_female_results", "median_results", 
                            "first_quartile_results", "threshold_results")) {
      c(0, 100)
    } else if (name %in% c("asymptote_male_results", "asymptote_female_results", "asymptote_results")) {
      c(0, 1) # Assuming asymptote values are between 0 and 1
    } else {
      range(param_data, na.rm = TRUE) # Default to data range
    }
    
    # Ensure xlim is finite and valid
    if (any(is.infinite(xlim))) {
      xlim <- c(min(param_data, na.rm = TRUE), max(param_data, na.rm = TRUE))
    }
    
    # Create the histogram
    hist(param_data,
         main = paste("Density Plot of", mod_name),
         xlab = mod_name,
         freq = FALSE,
         xlim = xlim,
         breaks = "Sturges", # Default break algorithm
         col = "lightblue" # Optional: color for the histogram
    )
  }
  
  # Reset to default single plot setting after plotting
  par(mfrow = c(1, 1))
}

plot_trace <- function(results, n_chains) {
  # Set up a grid for the plots based on the number of chains
  if (n_chains <= 3) {
    par(mfrow = c(n_chains * 2, 2)) # Up to 3 chains: 3 rows, 4 columns
  } else {
    par(mfrow = c(ceiling(n_chains), 4)) # More than 3 chains: 2 rows, 8 columns
  }
  
  for (chain_id in 1:n_chains) {
    if (!is.null(results[[chain_id]]$median_male_samples) || !is.null(results[[chain_id]]$median_female_samples)) {
      # Plot sex-specific parameters if available
      median_results <- results[[chain_id]]$median_male_samples
      threshold_results <- results[[chain_id]]$threshold_male_samples
      first_quartile_results <- results[[chain_id]]$first_quartile_male_samples
      asymptote_results <- results[[chain_id]]$asymptote_male_samples
      
      # Plot median, threshold, first quartile, and asymptote for males
      if (length(median_results) > 0) {
        plot(median_results, type = "l", main = paste("Chain", chain_id, "- Trace plot of Median - Male"), xlab = "Iteration", ylab = "Median")
      }
      if (length(threshold_results) > 0) {
        plot(threshold_results, type = "l", main = paste("Chain", chain_id, "- Trace plot of Threshold - Male"), xlab = "Iteration", ylab = "Threshold")
      }
      if (length(first_quartile_results) > 0) {
        plot(first_quartile_results, type = "l", main = paste("Chain", chain_id, "- Trace plot of First Quartile - Male"), xlab = "Iteration", ylab = "First Quartile")
      }
      if (length(asymptote_results) > 0) {
        plot(asymptote_results, type = "l", main = paste("Chain", chain_id, "- Trace plot of Asymptote - Male"), xlab = "Iteration", ylab = "Asymptote")
      }
      
      # Now plot for females
      median_results <- results[[chain_id]]$median_female_samples
      threshold_results <- results[[chain_id]]$threshold_female_samples
      first_quartile_results <- results[[chain_id]]$first_quartile_female_samples
      asymptote_results <- results[[chain_id]]$asymptote_female_samples
      
      if (length(median_results) > 0) {
        plot(median_results, type = "l", main = paste("Chain", chain_id, "- Trace plot of Median - Female"), xlab = "Iteration", ylab = "Median")
      }
      if (length(threshold_results) > 0) {
        plot(threshold_results, type = "l", main = paste("Chain", chain_id, "- Trace plot of Threshold - Female"), xlab = "Iteration", ylab = "Threshold")
      }
      if (length(first_quartile_results) > 0) {
        plot(first_quartile_results, type = "l", main = paste("Chain", chain_id, "- Trace plot of First Quartile - Female"), xlab = "Iteration", ylab = "First Quartile")
      }
      if (length(asymptote_results) > 0) {
        plot(asymptote_results, type = "l", main = paste("Chain", chain_id, "- Trace plot of Asymptote - Female"), xlab = "Iteration", ylab = "Asymptote")
      }
    } else {
      # Plot non-sex-specific parameters if sex-specific are not available
      median_results <- results[[chain_id]]$median_samples
      threshold_results <- results[[chain_id]]$threshold_samples
      first_quartile_results <- results[[chain_id]]$first_quartile_samples
      asymptote_results <- results[[chain_id]]$asymptote_samples
      
      if (length(median_results) > 0) {
        plot(median_results, type = "l", main = paste("Chain", chain_id, "- Trace plot of Median"), xlab = "Iteration", ylab = "Median")
      }
      if (length(threshold_results) > 0) {
        plot(threshold_results, type = "l", main = paste("Chain", chain_id, "- Trace plot of Threshold"), xlab = "Iteration", ylab = "Threshold")
      }
      if (length(first_quartile_results) > 0) {
        plot(first_quartile_results, type = "l", main = paste("Chain", chain_id, "- Trace plot of First Quartile"), xlab = "Iteration", ylab = "First Quartile")
      }
      if (length(asymptote_results) > 0) {
        plot(asymptote_results, type = "l", main = paste("Chain", chain_id, "- Trace plot of Asymptote"), xlab = "Iteration", ylab = "Asymptote")
      }
    }
  }
  
  # Reset the plotting layout
  par(mfrow = c(1, 1))
}

plot_traceSingle <- function(results) {
  par(mfrow = c(2, 2)) # Set up a grid for the plots
  
  if (!is.null(results$median_male_samples) || !is.null(results$median_female_samples)) {
    # Plot sex-specific parameters if available
    median_results <- results$median_male_samples
    threshold_results <- results$threshold_male_samples
    first_quartile_results <- results$first_quartile_male_samples
    asymptote_results <- results$asymptote_male_samples
    
    # Plot median, threshold, first quartile, and asymptote for males
    if (length(median_results) > 0) {
      plot(median_results, type = "l", main = "Trace plot of Median - Male", xlab = "Iteration", ylab = "Median")
    }
    if (length(threshold_results) > 0) {
      plot(threshold_results, type = "l", main = "Trace plot of Threshold - Male", xlab = "Iteration", ylab = "Threshold")
    }
    if (length(first_quartile_results) > 0) {
      plot(first_quartile_results, type = "l", main = "Trace plot of First Quartile - Male", xlab = "Iteration", ylab = "First Quartile")
    }
    if (length(asymptote_results) > 0) {
      plot(asymptote_results, type = "l", main = "Trace plot of Asymptote - Male", xlab = "Iteration", ylab = "Asymptote")
    }
    
    # Now plot for females
    median_results <- results$median_female_samples
    threshold_results <- results$threshold_female_samples
    first_quartile_results <- results$first_quartile_female_samples
    asymptote_results <- results$asymptote_female_samples
    
    if (length(median_results) > 0) {
      plot(median_results, type = "l", main = "Trace plot of Median - Female", xlab = "Iteration", ylab = "Median")
    }
    if (length(threshold_results) > 0) {
      plot(threshold_results, type = "l", main = "Trace plot of Threshold - Female", xlab = "Iteration", ylab = "Threshold")
    }
    if (length(first_quartile_results) > 0) {
      plot(first_quartile_results, type = "l", main = "Trace plot of First Quartile - Female", xlab = "Iteration", ylab = "First Quartile")
    }
    if (length(asymptote_results) > 0) {
      plot(asymptote_results, type = "l", main = "Trace plot of Asymptote - Female", xlab = "Iteration", ylab = "Asymptote")
    }
  } else {
    # Plot non-sex-specific parameters if sex-specific are not available
    median_results <- results$median_samples
    threshold_results <- results$threshold_samples
    first_quartile_results <- results$first_quartile_samples
    asymptote_results <- results$asymptote_samples
    
    if (length(median_results) > 0) {
      plot(median_results, type = "l", main = "Trace plot of Median", xlab = "Iteration", ylab = "Median")
    }
    if (length(threshold_results) > 0) {
      plot(threshold_results, type = "l", main = "Trace plot of Threshold", xlab = "Iteration", ylab = "Threshold")
    }
    if (length(first_quartile_results) > 0) {
      plot(first_quartile_results, type = "l", main = "Trace plot of First Quartile", xlab = "Iteration", ylab = "First Quartile")
    }
    if (length(asymptote_results) > 0) {
      plot(asymptote_results, type = "l", main = "Trace plot of Asymptote", xlab = "Iteration", ylab = "Asymptote")
    }
  }
  
  # Reset the plotting layout
  par(mfrow = c(1, 1))
}

plot_penetrance <- function(data, prob, max_age, sex = "NA") {
  if (prob <= 0 || prob >= 1) {
    stop("prob must be between 0 and 1")
  }
  
  # Check if sex-specific parameters are present
  if (!is.null(data$median_male_results) || !is.null(data$median_female_results)) {
    # Use sex-specific parameters
    params_male <- calculate_weibull_parameters(
      data$median_male_results,
      data$first_quartile_male_results,
      data$threshold_male_results
    )
    
    params_female <- calculate_weibull_parameters(
      data$median_female_results,
      data$first_quartile_female_results,
      data$threshold_female_results
    )
    
    alphas_male <- params_male$alpha
    betas_male <- params_male$beta
    thresholds_male <- data$threshold_male_results
    alphas_female <- params_female$alpha
    betas_female <- params_female$beta
    thresholds_female <- data$threshold_female_results
    
    asymptotes_male <- data$asymptote_male_results
    asymptotes_female <- data$asymptote_female_results
  } else {
    # Use non-sex-specific parameters
    params <- calculate_weibull_parameters(
      data$median_samples,
      data$first_quartile_samples,
      data$threshold_samples
    )
    
    alphas <- params$alpha
    betas <- params$beta
    thresholds <- data$threshold_samples
    asymptotes <- data$asymptote_samples
  }
  
  x_values <- seq(0, max_age, length.out = max_age + 1)
  
  plot_distribution <- function(alphas, betas, thresholds, asymptotes, x_values, prob, color, add = FALSE) {
    distributions <- mapply(function(alpha, beta, threshold, asymptote) {
      pweibull(x_values - threshold, shape = alpha, scale = beta) * asymptote
    }, alphas, betas, thresholds, asymptotes, SIMPLIFY = FALSE)
    
    distributions_matrix <- matrix(unlist(distributions), nrow = length(x_values), byrow = FALSE)
    mean_density <- rowMeans(distributions_matrix, na.rm = TRUE)
    ci_lower <- apply(distributions_matrix, 1, quantile, probs = (1 - prob) / 2, na.rm = TRUE)
    ci_upper <- apply(distributions_matrix, 1, quantile, probs = 1 - (1 - prob) / 2, na.rm = TRUE)
    
    if (!add) {
      plot(x_values, mean_density,
           type = "l", col = color,
           ylim = c(min(ci_lower, na.rm = TRUE), max(ci_upper, na.rm = TRUE)),
           xlab = "Age", ylab = "Cumulative Penetrance", main = "Penetrance Curve with Credible Interval - Cumulative Probability"
      )
    } else {
      lines(x_values, mean_density, col = color)
    }
    lines(x_values, ci_lower, col = color, lty = 2)
    lines(x_values, ci_upper, col = color, lty = 2)
    polygon(c(x_values, rev(x_values)), c(ci_lower, rev(ci_upper)), col = rgb(0, 0, 1, 0.1), border = NA)
  }
  
  if (!is.null(data$median_male_results) || !is.null(data$median_female_results)) {
    if (sex == "Male") {
      plot_distribution(alphas_male, betas_male, thresholds_male, asymptotes_male, x_values, prob, "blue")
      legend_text <- "Male"
    } else if (sex == "Female") {
      plot_distribution(alphas_female, betas_female, thresholds_female, asymptotes_female, x_values, prob, "red")
      legend_text <- "Female"
    } else {
      plot_distribution(alphas_male, betas_male, thresholds_male, asymptotes_male, x_values, prob, "blue")
      plot_distribution(alphas_female, betas_female, thresholds_female, asymptotes_female, x_values, prob, "red", add = TRUE)
      legend_text <- c("Male", "Female")
    }
  } else {
    # Plot for non-sex-specific parameters
    plot_distribution(alphas, betas, thresholds, asymptotes, x_values, prob, "green")
    legend_text <- "Overall"
  }
  
  legend("topleft",
         legend = legend_text,
         col = if (sex == "NA") c("blue", "red", "green") else c("blue", "red"),
         lty = c(1, 1),
         cex = 0.8
  )
}

plot_pdf <- function(data, prob, max_age, sex = "NA") {
  if (prob <= 0 || prob >= 1) {
    stop("prob must be between 0 and 1")
  }
  
  # Check if sex-specific parameters are present
  if (!is.null(data$median_male_results) || !is.null(data$median_female_results)) {
    # Use sex-specific parameters
    params_male <- calculate_weibull_parameters(
      data$median_male_results,
      data$first_quartile_male_results,
      data$threshold_male_results
    )
    
    params_female <- calculate_weibull_parameters(
      data$median_female_results,
      data$first_quartile_female_results,
      data$threshold_female_results
    )
    
    alphas_male <- params_male$alpha
    betas_male <- params_male$beta
    thresholds_male <- data$threshold_male_results
    alphas_female <- params_female$alpha
    betas_female <- params_female$beta
    thresholds_female <- data$threshold_female_results
    
    asymptotes_male <- data$asymptote_male_results
    asymptotes_female <- data$asymptote_female_results
  } else {
    # Use non-sex-specific parameters
    params <- calculate_weibull_parameters(
      data$median_samples,
      data$first_quartile_samples,
      data$threshold_samples
    )
    
    alphas <- params$alpha
    betas <- params$beta
    thresholds <- data$threshold_samples
    asymptotes <- data$asymptote_samples
  }
  
  x_values <- seq(0, max_age, length.out = max_age + 1)
  
  plot_pdf_distribution <- function(alphas, betas, thresholds, asymptotes, x_values, prob, color, add = FALSE) {
    pdf_distributions <- mapply(function(alpha, beta, threshold, asymptote) {
      dweibull(x_values - threshold, shape = alpha, scale = beta) * asymptote
    }, alphas, betas, thresholds, asymptotes, SIMPLIFY = FALSE)
    
    pdf_matrix <- matrix(unlist(pdf_distributions), nrow = length(x_values), byrow = FALSE)
    mean_density <- rowMeans(pdf_matrix, na.rm = TRUE)
    ci_lower <- apply(pdf_matrix, 1, quantile, probs = (1 - prob) / 2, na.rm = TRUE)
    ci_upper <- apply(pdf_matrix, 1, quantile, probs = 1 - (1 - prob) / 2, na.rm = TRUE)
    
    if (!add) {
      plot(x_values, mean_density,
           type = "l", col = color,
           ylim = c(min(ci_lower, na.rm = TRUE), max(ci_upper, na.rm = TRUE)),
           xlab = "Age", ylab = "Probability Density", main = "Penetrance Curve with Credible Interval - Probability Distribution"
      )
    } else {
      lines(x_values, mean_density, col = color)
    }
    lines(x_values, ci_lower, col = color, lty = 2)
    lines(x_values, ci_upper, col = color, lty = 2)
    polygon(c(x_values, rev(x_values)), c(ci_lower, rev(ci_upper)), col = rgb(0, 0, 1, 0.1), border = NA)
  }
  
  if (!is.null(data$median_male_results) || !is.null(data$median_female_results)) {
    # Plot for sex-specific parameters
    if (sex == "Male") {
      plot_pdf_distribution(alphas_male, betas_male, thresholds_male, asymptotes_male, x_values, prob, "blue")
      legend_text <- "Male"
    } else if (sex == "Female") {
      plot_pdf_distribution(alphas_female, betas_female, thresholds_female, asymptotes_female, x_values, prob, "red")
      legend_text <- "Female"
    } else {
      plot_pdf_distribution(alphas_male, betas_male, thresholds_male, asymptotes_male, x_values, prob, "blue")
      plot_pdf_distribution(alphas_female, betas_female, thresholds_female, asymptotes_female, x_values, prob, "red", add = TRUE)
      legend_text <- c("Male", "Female")
    }
  } else {
    # Plot for non-sex-specific parameters
    plot_pdf_distribution(alphas, betas, thresholds, asymptotes, x_values, prob, "green")
    legend_text <- "Overall"
  }
  
  legend("topleft",
         legend = legend_text,
         col = if (sex == "NA") c("blue", "red", "green") else c("blue", "red"),
         lty = c(1, 1),
         cex = 0.8
  )
}

combine_chains_noSex <- function(results) {
  list(
    median_results = do.call(c, lapply(results, function(x) x$median_samples)),
    threshold_results = do.call(c, lapply(results, function(x) x$threshold_samples)),
    first_quartile_results = do.call(c, lapply(results, function(x) x$first_quartile_samples)),
    asymptote_results = do.call(c, lapply(results, function(x) x$asymptote_samples)),
    loglikelihood_current_results = do.call(c, lapply(results, function(x) x$loglikelihood_current)),
    loglikelihood_proposal_results = do.call(c, lapply(results, function(x) x$loglikelihood_proposal)),
    acceptance_ratio_results = do.call(c, lapply(results, function(x) x$acceptance_ratio)),
    median_proposals = do.call(c, lapply(results, function(x) x$median_proposals)),
    threshold_proposals = do.call(c, lapply(results, function(x) x$threshold_proposals)),
    first_quartile_proposals = do.call(c, lapply(results, function(x) x$first_quartile_proposals)),
    asymptote_proposals = do.call(c, lapply(results, function(x) x$asymptote_proposals))
  )
}

#' Generate Summary
#' @description Function to generate summary statistics
#'
#' @param data A list with combined results.
#'
#' @return A summary data frame containing Median, threshold, First Quartile, and Asymptote Value.
#' @export
generate_summary_noSex <- function(data) {
  summary_data <- data.frame(
    Median = data$median_results,
    Threshold = data$threshold_results,
    First_Quartile = data$first_quartile_results,
    Asymptote = data$asymptote_results
  )
  print(summary(summary_data))
  return(summary(summary_data))
}
