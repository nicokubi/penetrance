#' Combine Chains
#' Function to combine the posterior samples from the multiple chains.
#'
#' @param results A list of MCMC chain results.
#'
#' @return A list with combined results, including median, threshold, first quartile, and asymptote values.
#' 
combine_chains <- function(results) {
  list(
    median_male_results = do.call(c, lapply(results, function(x) x$median_male_samples)),
    median_female_results = do.call(c, lapply(results, function(x) x$median_female_samples)),
    threshold_male_results = do.call(c, lapply(results, function(x) x$threshold_male_samples)),
    threshold_female_results = do.call(c, lapply(results, function(x) x$threshold_female_samples)),
    first_quartile_male_results = do.call(c, lapply(results, function(x) x$first_quartile_male_samples)),
    first_quartile_female_results = do.call(c, lapply(results, function(x) x$first_quartile_female_samples)),
    asymptote_male_results = do.call(c, lapply(results, function(x) x$asymptote_male_samples)),
    asymptote_female_results = do.call(c, lapply(results, function(x) x$asymptote_female_samples)),
    loglikelihood_current_results = do.call(c, lapply(results, function(x) x$loglikelihood_current)),
    loglikelihood_proposal_results = do.call(c, lapply(results, function(x) x$loglikelihood_proposal)),
    log_acceptance_ratio_results = do.call(c, lapply(results, function(x) x$log_acceptance_ratio)),
    median_male_proposals = do.call(c, lapply(results, function(x) x$median_male_proposals)),
    median_female_proposals = do.call(c, lapply(results, function(x) x$median_female_proposals)),
    threshold_male_proposals = do.call(c, lapply(results, function(x) x$threshold_male_proposals)),
    threshold_female_proposals = do.call(c, lapply(results, function(x) x$threshold_female_proposals)),
    first_quartile_male_proposals = do.call(c, lapply(results, function(x) x$first_quartile_male_proposals)),
    first_quartile_female_proposals = do.call(c, lapply(results, function(x) x$first_quartile_female_proposals)),
    asymptote_male_proposals = do.call(c, lapply(results, function(x) x$asymptote_male_proposals)),
    asymptote_female_proposals = do.call(c, lapply(results, function(x) x$asymptote_female_proposals))
  )
}

#' Generate Summary
#' @description Function to generate summary statistics
#'
#' @param data A list with combined results.
#'
#' @return A summary data frame containing Median, threshold, First Quartile, and Asymptote Value.
#' @export
generate_summary <- function(data) {
  summary_data <- data.frame(
    Median_Male = data$median_male_results,
    Median_Female = data$median_female_results,
    Threshold_Male = data$threshold_male_results,
    Threshold_Female = data$threshold_female_results,
    First_Quartile_Male = data$first_quartile_male_results,
    First_Quartile_Female = data$first_quartile_female_results,
    Asymptote_Male = data$asymptote_male_results,
    Asymptote_Female = data$asymptote_female_results
  )
  print(summary(summary_data))
  return(summary(summary_data))
}

#' Generate Posterior Density Plots
#' 
#' Generates histograms of the posterior samples for the different parameters
#'
#' @param data A list with combined results.
#' 
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

#' Plot Trace
#' @param results A list of MCMC chain results.
#' @param n_chains The number of chains.
#' @export
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

# Running mean calculation
running_mean <- function(res) {
  n <- length(res)
  running_mean <- numeric(n)
  for (i in 1:n) {
    running_mean[i] <- mean(res[1:i])
  }
  return(running_mean)
}

# Running var calculation
running_variance <- function(res) {
  n <- length(res)
  running_var <- numeric(n)
  for (i in 1:n) {
    running_var[i] <- var(res[1:i])
  }
  return(running_var)
}

#' Print Rejection Rates
#' @param results A list of MCMC chain results.
#' @export
printRejectionRates <- function(results) {
  rejection_rates <- sapply(results, function(x) x$rejection_rate)
  cat("Rejection rates: ", rejection_rates, "\n")
}

#' Apply Burn-In
#'
#' @param results A list of MCMC chain results.
#' @param burn_in The fraction roportion of results to discard as burn-in (0 to 1). The default is no burn-in, burn_in=0.
#'
#' @return A list of results with burn-in applied.
#'
apply_burn_in <- function(results, burn_in) {
  # Ensure 'results' is a list and has at least one chain
  if (!is.list(results) || length(results) < 1) {
    stop("results must be a list with at least one chain.")
  }

  # Ensure 'burn_in' is numeric and between 0 and 1
  if (!is.numeric(burn_in) || burn_in <= 0 || burn_in >= 1) {
    stop("burn_in must be a numeric value between 0 and 1.")
  }

  # Function to perform burn-in on a single chain (list of numeric vectors)
  burn_in_chain <- function(chain, burn_in) {
    lapply(chain, function(param_results) {
      n_results <- length(param_results)
      burn_in_count <- round(n_results * burn_in)
      if (burn_in_count >= n_results) {
        stop("burn_in_count must be less than the total number of results (n_results).")
      }
      param_results[(burn_in_count + 1):n_results]
    })
  }

  # Apply burn-in to all results
  lapply(results, function(chain) {
    burn_in_chain(chain, burn_in)
  })
}

##' Apply Thinning
#'
#' @param results A list of MCMC chain results.
#' @param thinning_factor The factor by which to thin the results (positive integer). The default thinning factor is 1, which implies no thinning.
#'
#' @return A list of results with thinning applied.
#'
apply_thinning <- function(results, thinning_factor) {
  if (!is.numeric(thinning_factor) || thinning_factor <= 0 || thinning_factor != round(thinning_factor)) {
    stop("thinning_factor must be a positive integer.")
  }
  thinning_factor <- as.integer(thinning_factor)
  
  # Define a function to thin each element of the chain
  thin_list <- function(chain, factor) {
    lapply(chain, function(param) {
      if (is.vector(param)) {
        # If param is a vector, select every nth element
        return(param[seq(1, length(param), by = factor)])
      } else if (is.matrix(param)) {
        # If param is a matrix, select every nth row
        return(param[seq(1, nrow(param), by = factor), , drop = FALSE])
      } else if (is.list(param)) {
        # If param is a list, recursively thin each element of the list
        return(lapply(param, function(sub_param) {
          if (is.vector(sub_param)) {
            return(sub_param[seq(1, length(sub_param), by = factor)])
          } else if (is.matrix(sub_param)) {
            return(sub_param[seq(1, nrow(sub_param), by = factor), , drop = FALSE])
          } else {
            return(sub_param) # If it's not a vector/matrix, return as is
          }
        }))
      } else {
        return(param) # If it's not a vector/matrix/list, return as is
      }
    })
  }
  
  # Apply thinning to the results
  thinned_results <- lapply(results, function(chain) {
    thin_list(chain, thinning_factor)
  })
  
  return(thinned_results)
}

# Example usage:
# thinned_results <- apply_thinning(out_OC_PALB2$results, 5)


#' Plot Weibull Distribution with Credible Intervals
#'
#' This function plots the Weibull distribution with credible intervals for the given data. 
#' It allows for visualization of penetrance curves for individuals based on their genetic 
#' and demographic information.
#'
#' @param data Data frame, containing individual demographic and genetic information. Must include columns for 'sex', 'age', 'aff' (affection status), and 'geno' (genotype).
#' @param prob Numeric, the probability level for the credible intervals. Must be between 0 and 1.
#' @param max_age Integer, the maximum age considered in the analysis.
#' @param sex Character, specifying the sex of the individuals for the plot ("Male", "Female", or "NA" for not applicable). Default is "NA".
#'
#' @return A plot showing the Weibull distribution with credible intervals.
#' @export
plot_penetrance <- function(data, prob, max_age, sex = "NA") {
  if (prob <= 0 || prob >= 1) {
    stop("prob must be between 0 and 1")
  }
  
  # Check if sex-specific parameters are present
  sex_specific <- !is.null(data$median_male_results) && !is.null(data$median_female_results)
  
  if (sex_specific) {
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
      data$median_results,
      data$first_quartile_results,
      data$threshold_results
    )
    
    alphas <- params$alpha
    betas <- params$beta
    thresholds <- data$threshold_results
    asymptotes <- data$asymptote_results
  }
  
  x_values <- seq(0, max_age, length.out = max_age + 1)
  
  calculate_ylim <- function(alphas, betas, thresholds, asymptotes, x_values, prob) {
    distributions <- mapply(function(alpha, beta, threshold, asymptote) {
      pweibull(x_values - threshold, shape = alpha, scale = beta) * asymptote
    }, alphas, betas, thresholds, asymptotes, SIMPLIFY = FALSE)
    
    distributions_matrix <- matrix(unlist(distributions), nrow = length(x_values), byrow = FALSE)
    ci_lower <- apply(distributions_matrix, 1, quantile, probs = (1 - prob) / 2, na.rm = TRUE)
    ci_upper <- apply(distributions_matrix, 1, quantile, probs = 1 - (1 - prob) / 2, na.rm = TRUE)
    
    return(c(min(ci_lower, na.rm = TRUE), max(ci_upper, na.rm = TRUE)))
  }
  
  plot_distribution <- function(alphas, betas, thresholds, asymptotes, x_values, prob, color, add = FALSE, ylim = NULL) {
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
           ylim = ylim,
           xlab = "Age", ylab = "Cumulative Penetrance", 
           main = "Penetrance Curve with Credible Interval - Cumulative Probability"
      )
    } else {
      lines(x_values, mean_density, col = color)
    }
    lines(x_values, ci_lower, col = color, lty = 2)
    lines(x_values, ci_upper, col = color, lty = 2)
    polygon(c(x_values, rev(x_values)), c(ci_lower, rev(ci_upper)), col = adjustcolor(color, alpha.f = 0.1), border = NA)
  }
  
  if (sex_specific) {
    # Calculate y-limits for both male and female distributions without plotting
    ylim_male <- calculate_ylim(alphas_male, betas_male, thresholds_male, asymptotes_male, x_values, prob)
    ylim_female <- calculate_ylim(alphas_female, betas_female, thresholds_female, asymptotes_female, x_values, prob)
    
    # Combine y-limits
    combined_ylim <- c(min(ylim_male[1], ylim_female[1]), max(ylim_male[2], ylim_female[2]))
    
    # Plot for sex-specific parameters with combined y-limits
    if (sex == "Male") {
      plot_distribution(alphas_male, betas_male, thresholds_male, asymptotes_male, x_values, prob, "blue", add = FALSE, ylim = combined_ylim)
      legend_text <- "Male"
    } else if (sex == "Female") {
      plot_distribution(alphas_female, betas_female, thresholds_female, asymptotes_female, x_values, prob, "red", add = FALSE, ylim = combined_ylim)
      legend_text <- "Female"
    } else {
      plot_distribution(alphas_male, betas_male, thresholds_male, asymptotes_male, x_values, prob, "blue", add = FALSE, ylim = combined_ylim)
      plot_distribution(alphas_female, betas_female, thresholds_female, asymptotes_female, x_values, prob, "red", add = TRUE, ylim = combined_ylim)
      legend_text <- c("Male", "Female")
    }
  } else {
    # Plot for non-sex-specific parameters
    plot_distribution(alphas, betas, thresholds, asymptotes, x_values, prob, "green", add = FALSE)
    legend_text <- "Overall"
  }
  
  legend("topleft",
         legend = legend_text,
         col = if (sex_specific && sex == "NA") c("blue", "red") else "green",
         lty = c(1, 1),
         cex = 0.8
  )
}

#' Plot Weibull Probability Density Function with Credible Intervals
#'
#' This function plots the Weibull PDF with credible intervals for the given data. 
#' It allows for visualization of density curves for individuals based on their genetic 
#' and demographic information.
#'
#' @param data Data frame, containing individual demographic and genetic information. Must include columns for 'sex', 'age', 'aff' (affection status), and 'geno' (genotype).
#' @param prob Numeric, the probability level for the credible intervals. Must be between 0 and 1.
#' @param max_age Integer, the maximum age considered in the analysis.
#' @param sex Character, specifying the sex of the individuals for the plot ("Male", "Female", or "NA" for not applicable). Default is "NA".
#'
#' @return A plot showing the Weibull PDF with credible intervals.
#' @export
plot_pdf <- function(data, prob, max_age, sex = "NA") {
  if (prob <= 0 || prob >= 1) {
    stop("prob must be between 0 and 1")
  }
  
  # Check if sex-specific parameters are present
  sex_specific <- !is.null(data$median_male_results) && !is.null(data$median_female_results)
  
  if (sex_specific) {
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
      data$median_results,
      data$first_quartile_results,
      data$threshold_results
    )
    
    alphas <- params$alpha
    betas <- params$beta
    thresholds <- data$threshold_results
    asymptotes <- data$asymptote_results
  }
  
  x_values <- seq(0, max_age, length.out = max_age + 1)
  
  calculate_ylim <- function(alphas, betas, thresholds, asymptotes, x_values, prob) {
    pdf_distributions <- mapply(function(alpha, beta, threshold, asymptote) {
      dweibull(x_values - threshold, shape = alpha, scale = beta) * asymptote
    }, alphas, betas, thresholds, asymptotes, SIMPLIFY = FALSE)
    
    pdf_matrix <- matrix(unlist(pdf_distributions), nrow = length(x_values), byrow = FALSE)
    ci_lower <- apply(pdf_matrix, 1, quantile, probs = (1 - prob) / 2, na.rm = TRUE)
    ci_upper <- apply(pdf_matrix, 1, quantile, probs = 1 - (1 - prob) / 2, na.rm = TRUE)
    
    return(c(min(ci_lower, na.rm = TRUE), max(ci_upper, na.rm = TRUE)))
  }
  
  plot_pdf_distribution <- function(alphas, betas, thresholds, asymptotes, x_values, prob, color, add = FALSE, ylim = NULL) {
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
           ylim = ylim,
           xlab = "Age", ylab = "Probability Density", 
           main = "Penetrance Curve with Credible Interval - Probability Distribution"
      )
    } else {
      lines(x_values, mean_density, col = color)
    }
    lines(x_values, ci_lower, col = color, lty = 2)
    lines(x_values, ci_upper, col = color, lty = 2)
    polygon(c(x_values, rev(x_values)), c(ci_lower, rev(ci_upper)), col = adjustcolor(color, alpha.f = 0.1), border = NA)
  }
  
  if (sex_specific) {
    # Calculate y-limits for both male and female distributions without plotting
    ylim_male <- calculate_ylim(alphas_male, betas_male, thresholds_male, asymptotes_male, x_values, prob)
    ylim_female <- calculate_ylim(alphas_female, betas_female, thresholds_female, asymptotes_female, x_values, prob)
    
    # Combine y-limits
    combined_ylim <- c(min(ylim_male[1], ylim_female[1]), max(ylim_male[2], ylim_female[2]))
    
    # Plot for sex-specific parameters with combined y-limits
    if (sex == "Male") {
      plot_pdf_distribution(alphas_male, betas_male, thresholds_male, asymptotes_male, x_values, prob, "blue", add = FALSE, ylim = combined_ylim)
      legend_text <- "Male"
    } else if (sex == "Female") {
      plot_pdf_distribution(alphas_female, betas_female, thresholds_female, asymptotes_female, x_values, prob, "red", add = FALSE, ylim = combined_ylim)
      legend_text <- "Female"
    } else {
      plot_pdf_distribution(alphas_male, betas_male, thresholds_male, asymptotes_male, x_values, prob, "blue", add = FALSE, ylim = combined_ylim)
      plot_pdf_distribution(alphas_female, betas_female, thresholds_female, asymptotes_female, x_values, prob, "red", add = TRUE, ylim = combined_ylim)
      legend_text <- c("Male", "Female")
    }
  } else {
    # Plot for non-sex-specific parameters
    plot_pdf_distribution(alphas, betas, thresholds, asymptotes, x_values, prob, "green", add = FALSE)
    legend_text <- "Overall"
  }
  
  legend("topleft",
         legend = legend_text,
         col = if (sex_specific && sex == "NA") c("blue", "red") else "green",
         lty = c(1, 1),
         cex = 0.8
  )
}

#' Combine Chains for Non-Sex-Specific Estimation
#'
#' Combines the posterior samples from multiple MCMC chains for non-sex-specific estimations.
#'
#' @param results A list of MCMC chain results, where each element contains posterior samples of parameters.
#'
#' @return A list with combined results, including samples for median, threshold, first quartile, asymptote values, 
#' log-likelihoods, and log-acceptance ratios.
#' 
#' @export
combine_chains_noSex <- function(results) {
  list(
    median_results = do.call(c, lapply(results, function(x) x$median_samples)),
    threshold_results = do.call(c, lapply(results, function(x) x$threshold_samples)),
    first_quartile_results = do.call(c, lapply(results, function(x) x$first_quartile_samples)),
    asymptote_results = do.call(c, lapply(results, function(x) x$asymptote_samples)),
    loglikelihood_current_results = do.call(c, lapply(results, function(x) x$loglikelihood_current)),
    loglikelihood_proposal_results = do.call(c, lapply(results, function(x) x$loglikelihood_proposal)),
    log_acceptance_ratio_results = do.call(c, lapply(results, function(x) x$log_acceptance_ratio)),
    median_proposals = do.call(c, lapply(results, function(x) x$median_proposals)),
    threshold_proposals = do.call(c, lapply(results, function(x) x$threshold_proposals)),
    first_quartile_proposals = do.call(c, lapply(results, function(x) x$first_quartile_proposals)),
    asymptote_proposals = do.call(c, lapply(results, function(x) x$asymptote_proposals))
  )
}

#' Generate Summary for Non-Sex-Specific Estimation
#'
#' Generates summary statistics for the combined MCMC results for non-sex-specific estimations.
#'
#' @param data A list containing combined results of MCMC chains, typically the output of `combine_chains_noSex`.
#'
#' @return A summary data frame containing median, threshold, first quartile, and asymptote values.
#' 
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

#' Plot Autocorrelation for Multiple MCMC Chains (Posterior Samples)
#'
#' This function plots the autocorrelation for sex-specific or non-sex-specific posterior samples across multiple MCMC chains. 
#' It defaults to key parameters like `asymptote_male_samples`, `asymptote_female_samples`, etc.
#'
#' @param results A list of MCMC chain results.
#' @param n_chains The number of chains.
#' @param max_lag Integer, the maximum lag to be considered for the autocorrelation plot. Default is 50.
#'
#' @return A series of autocorrelation plots for each chain.
#' @export
plot_acf <- function(results, n_chains, max_lag = 50) {
  # Set up a grid for the plots based on the number of chains
  if (n_chains <= 3) {
    par(mfrow = c(n_chains * 2, 2))  # Up to 3 chains: 3 rows, 4 columns
  } else {
    par(mfrow = c(ceiling(n_chains), 4))  # More than 3 chains: grid layout
  }
  
  # Loop through each chain
  for (chain_id in 1:n_chains) {
    if (!is.null(results[[chain_id]]$median_male_samples) || !is.null(results[[chain_id]]$median_female_samples)) {
      # Plot ACF for sex-specific parameters if available
      median_results <- results[[chain_id]]$median_male_samples
      threshold_results <- results[[chain_id]]$threshold_male_samples
      first_quartile_results <- results[[chain_id]]$first_quartile_male_samples
      asymptote_results <- results[[chain_id]]$asymptote_male_samples
      
      # ACF plot for male parameters
      if (length(median_results) > 0) {
        acf(median_results, lag.max = max_lag, main = paste("Chain", chain_id, "- ACF of Median - Male"))
      }
      if (length(threshold_results) > 0) {
        acf(threshold_results, lag.max = max_lag, main = paste("Chain", chain_id, "- ACF of Threshold - Male"))
      }
      if (length(first_quartile_results) > 0) {
        acf(first_quartile_results, lag.max = max_lag, main = paste("Chain", chain_id, "- ACF of First Quartile - Male"))
      }
      if (length(asymptote_results) > 0) {
        acf(asymptote_results, lag.max = max_lag, main = paste("Chain", chain_id, "- ACF of Asymptote - Male"))
      }
      
      # Now plot for female parameters
      median_results <- results[[chain_id]]$median_female_samples
      threshold_results <- results[[chain_id]]$threshold_female_samples
      first_quartile_results <- results[[chain_id]]$first_quartile_female_samples
      asymptote_results <- results[[chain_id]]$asymptote_female_samples
      
      if (length(median_results) > 0) {
        acf(median_results, lag.max = max_lag, main = paste("Chain", chain_id, "- ACF of Median - Female"))
      }
      if (length(threshold_results) > 0) {
        acf(threshold_results, lag.max = max_lag, main = paste("Chain", chain_id, "- ACF of Threshold - Female"))
      }
      if (length(first_quartile_results) > 0) {
        acf(first_quartile_results, lag.max = max_lag, main = paste("Chain", chain_id, "- ACF of First Quartile - Female"))
      }
      if (length(asymptote_results) > 0) {
        acf(asymptote_results, lag.max = max_lag, main = paste("Chain", chain_id, "- ACF of Asymptote - Female"))
      }
    } else {
      # Plot ACF for non-sex-specific parameters if sex-specific are not available
      median_results <- results[[chain_id]]$median_samples
      threshold_results <- results[[chain_id]]$threshold_samples
      first_quartile_results <- results[[chain_id]]$first_quartile_samples
      asymptote_results <- results[[chain_id]]$asymptote_samples
      
      if (length(median_results) > 0) {
        acf(median_results, lag.max = max_lag, main = paste("Chain", chain_id, "- ACF of Median"))
      }
      if (length(threshold_results) > 0) {
        acf(threshold_results, lag.max = max_lag, main = paste("Chain", chain_id, "- ACF of Threshold"))
      }
      if (length(first_quartile_results) > 0) {
        acf(first_quartile_results, lag.max = max_lag, main = paste("Chain", chain_id, "- ACF of First Quartile"))
      }
      if (length(asymptote_results) > 0) {
        acf(asymptote_results, lag.max = max_lag, main = paste("Chain", chain_id, "- ACF of Asymptote"))
      }
    }
  }
  
  # Reset the plotting layout
  par(mfrow = c(1, 1))
}

#' Plot Log-Likelihood for Multiple MCMC Chains
#'
#' This function plots the log-likelihood values across iterations for multiple MCMC chains. 
#' It helps visualize the convergence of the chains based on the log-likelihood values.
#'
#' @param results A list of MCMC chain results, each containing the `loglikelihood_current` values.
#' @param n_chains The number of chains.
#'
#' @return A series of log-likelihood plots for each chain.
#' @export
plot_loglikelihood <- function(results, n_chains) {
  # Set up a grid for the plots based on the number of chains
  if (n_chains <= 3) {
    par(mfrow = c(n_chains, 1))  # Up to 3 chains: stacked vertically
  } else {
    par(mfrow = c(ceiling(n_chains / 2), 2))  # More than 3 chains: grid layout
  }
  
  # Loop through each chain
  for (chain_id in 1:n_chains) {
    # Check if the loglikelihood_current exists in the chain results
    if (is.null(results[[chain_id]]$loglikelihood_current)) {
      stop(paste("loglikelihood_current not found in chain", chain_id))
    }
    
    # Extract the log-likelihood values
    loglikelihood_values <- results[[chain_id]]$loglikelihood_current
    
    # Plot the log-likelihood values
    plot(loglikelihood_values, type = "l", col = "blue",
         main = paste("Chain", chain_id, "- Log-Likelihood"),
         xlab = "Iteration", ylab = "Log-Likelihood",
         ylim = range(loglikelihood_values, na.rm = TRUE))
    
    # Add a grid for better readability
    grid()
  }
  
  # Reset the plotting layout
  par(mfrow = c(1, 1))
}