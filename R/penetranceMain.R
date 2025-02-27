#' Bayesian Inference using Independent Metropolis-Hastings for Penetrance Estimation
#'
#' This function implements the Independent Metropolis-Hastings algorithm for Bayesian
#' penetrance estimation of cancer risk. It utilizes parallel computing to run multiple
#' chains and provides various options for analyzing and visualizing the results.
#'
#' @param pedigree A data frame containing the pedigree data in the required format. It should include the following columns:
#'   - `PedigreeID`: A numeric value representing the unique identifier for each family. There should be no duplicated entries.
#'   - `ID`: A numeric value representing the unique identifier for each individual. There should be no duplicated entries.
#'   - `Sex`: A numeric value where `0` indicates female and `1` indicates male. Unknown sex needs to be coded as `NA`. 
#'   - `MotherID`: A numeric value representing the unique identifier for an individual's mother.
#'   - `FatherID`: A numeric value representing the unique identifier for an individual's father.
#'   - `isProband`: A numeric value where `1` indicates the individual is a proband and `0` otherwise.
#'   - `CurAge`: A numeric value indicating the age of censoring (current age if the person is alive or age at death if the person is deceased). Allowed ages range from `1` to `94`. Unknown ages can be left empty or coded as `NA`. 
#'   - `isAff`: A numeric value indicating the affection status of cancer, with `1` for diagnosed individuals,
#'     `0` for unaffected individuals, and `NA` for unknown status.
#'   - `Age`: A numeric value indicating the age of cancer diagnosis, encoded as `NA` if the individual was not diagnosed. Allowed ages range from `1` to `94`. Unknown ages can be left empty or coded as `NA`. 
#'   - `geno`: A column for germline testing or tumor marker testing results. Positive results should be coded as `1`, negative results as `0`, and unknown results as `NA` or left empty.
#' @param twins A list specifying identical twins or triplets in the family. For example, to indicate that "ora024" and "ora027" are identical twins, and "aey063" and "aey064" are identical twins, use the following format: `twins <- list(c("ora024", "ora027"), c("aey063", "aey064"))`.
#' @param n_chains Integer, the number of chains for parallel computation. Default is 1.
#' @param n_iter_per_chain Integer, the number of iterations for each chain. Default is 10000.
#' @param ncores Integer, the number of cores for parallel computation. Default is 6.
#' @param baseline_data Data for the baseline risk estimates (probability of developing cancer), such as population-level risk from a cancer registry. Default data, for exemplary purposes, is for Colorectal cancer from the SEER database.
#' @param max_age Integer, the maximum age considered for analysis. Default is 94.
#' @param remove_proband Logical, indicating whether to remove probands from the analysis. Default is FALSE.
#' @param age_imputation Logical, indicating whether to perform age imputation. Default is FALSE.
#' @param median_max Logical, indicating whether to use the baseline median age or `max_age` as an upper bound for the median proposal. Default is TRUE.
#' @param BaselineNC Logical, indicating that the non-carrier penetrance is assumed to be the baseline penetrance. Default is TRUE.
#' @param var Numeric vector, variances for the proposal distribution in the Metropolis-Hastings algorithm. Default is `c(0.1, 0.1, 2, 2, 5, 5, 5, 5)`.
#' @param burn_in Numeric, the fraction of results to discard as burn-in (0 to 1). Default is 0 (no burn-in).
#' @param thinning_factor Integer, the factor by which to thin the results. Default is 1 (no thinning).
#' @param imp_interval Integer, the interval at which age imputation should be performed when age_imputation = TRUE.
#' @param distribution_data Data for generating prior distributions.
#' @param prev Numeric, prevalence of the carrier status. Default is 0.0001.
#' @param sample_size Optional numeric, sample size for distribution generation.
#' @param ratio Optional numeric, ratio parameter for distribution generation.
#' @param prior_params List, parameters for prior distributions.
#' @param risk_proportion Numeric, proportion of risk for distribution generation.
#' @param summary_stats Logical, indicating whether to include summary statistics in the output. Default is TRUE.
#' @param rejection_rates Logical, indicating whether to include rejection rates in the output. Default is TRUE.
#' @param density_plots Logical, indicating whether to include density plots in the output. Default is TRUE.
#' @param plot_trace Logical, indicating whether to include trace plots in the output. Default is TRUE.
#' @param penetrance_plot Logical, indicating whether to include penetrance plots in the output. Default is TRUE.
#' @param penetrance_plot_pdf Logical, indicating whether to include PDF plots in the output. Default is TRUE.
#' @param plot_loglikelihood Logical, indicating whether to include log-likelihood plots in the output. Default is TRUE.
#' @param plot_acf Logical, indicating whether to include autocorrelation function (ACF) plots for posterior samples. Default is TRUE.
#' @param probCI Numeric, probability level for credible intervals in penetrance plots. Must be between 0 and 1. Default is 0.95.
#' @param sex_specific Logical, indicating whether to use sex-specific parameters in the analysis. Default is TRUE.
#'
#' @return A list containing combined results from all chains, including optional statistics and plots.
#'
#' @importFrom stats rbeta runif
#' @importFrom parallel makeCluster stopCluster parLapply
#' 
#' @examples
#' # Create example baseline data (simplified for demonstration)
#' baseline_data_default <- data.frame(
#'   Age = 1:94,
#'   Female = rep(0.01, 94),
#'   Male = rep(0.01, 94)
#' )
#'
#' # Create example distribution data
#' distribution_data_default <- data.frame(
#'   Age = 1:94,
#'   Risk = rep(0.01, 94)
#' )
#'
#' # Create example prior parameters
#' prior_params_default <- list(
#'   shape = 2,
#'   scale = 50
#' )
#'
#' # Create example risk proportion
#' risk_proportion_default <- 0.5
#'
#' # Create a simple example pedigree
#' example_pedigree <- data.frame(
#'   PedigreeID = rep(1, 4),
#'   ID = 1:4,
#'   Sex = c(1, 0, 1, 0),  # 1 for male, 0 for female
#'   MotherID = c(NA, NA, 2, 2),
#'   FatherID = c(NA, NA, 1, 1),
#'   isProband = c(0, 0, 1, 0),
#'   CurAge = c(70, 68, 45, 42),
#'   isAff = c(0, 0, 1, 0),
#'   Age = c(NA, NA, 40, NA),
#'   geno = c(NA, NA, 1, NA)
#' )
#' 
#' # Basic usage with minimal iterations
#' result <- penetrance(
#'   pedigree = list(example_pedigree),
#'   n_chains = 1,
#'   n_iter_per_chain = 10,  # Very small number for example
#'   ncores = 1,             # Single core for example
#'   summary_stats = TRUE,
#'   plot_trace = FALSE,     # Disable plots for quick example
#'   density_plots = FALSE,
#'   penetrance_plot = FALSE,
#'   penetrance_plot_pdf = FALSE,
#'   plot_loglikelihood = FALSE,
#'   plot_acf = FALSE
#' )
#' 
#' # View basic results
#' head(result$summary_stats)
#' 
#' @export
penetrance <- function(pedigree,
                       twins = NULL,
                       n_chains = 1,
                       n_iter_per_chain = 10000,
                       ncores = 6,
                       max_age = 94,
                       baseline_data = baseline_data_default,
                       remove_proband = FALSE,
                       age_imputation = FALSE,
                       median_max = TRUE,
                       BaselineNC = TRUE,
                       var = c(0.1, 0.1, 2, 2, 5, 5, 5, 5),
                       burn_in = 0,
                       thinning_factor = 1,
                       imp_interval = 100,
                       distribution_data = distribution_data_default,
                       prev = 0.0001,
                       sample_size = NULL,
                       ratio = NULL,
                       prior_params = prior_params_default,
                       risk_proportion = risk_proportion_default,
                       summary_stats = TRUE,
                       rejection_rates = TRUE,
                       density_plots = TRUE,
                       plot_trace = TRUE,
                       penetrance_plot = TRUE,
                       penetrance_plot_pdf = TRUE,
                       plot_loglikelihood = TRUE,
                       plot_acf = TRUE,
                       probCI = 0.95,
                       sex_specific = TRUE) {
  # Validate inputs
  if (missing(pedigree) || !is.list(pedigree) || length(pedigree) == 0) {
    stop("Error: 'pedigree' parameter is missing or invalid. Please provide a non-empty list of pedigrees.")
  }

  # Check pedigree data structure and content for each pedigree in the list
  required_columns <- c("PedigreeID", "ID", "Sex", "MotherID", "FatherID", "isProband", "CurAge", "isAff", "Age", "geno")
  for (i in seq_along(pedigree)) {
    if (!all(required_columns %in% colnames(pedigree[[i]]))) {
      stop(paste("Error: Pedigree", i, "is missing one or more required columns."))
    }

    # Check for NA values in critical columns
    critical_columns <- c("PedigreeID", "ID")
    for (col in critical_columns) {
      if (any(is.na(pedigree[[i]][[col]]))) {
        stop(paste("Error: NA values found in the", col, "column of pedigree", i, ". Please check your data."))
      }
    }

    # Check isAff values
    if (!all(pedigree[[i]]$isAff %in% c(0, 1, NA))) {
      stop(paste("Error: 'isAff' column in pedigree", i, "should only contain 0, 1, or NA."))
    }
  }
  if (missing(n_chains) || !is.numeric(n_chains) || n_chains <= 0) {
    stop("Error: 'n_chains' parameter is missing or invalid. Please specify a positive integer.")
  }
  if (missing(n_iter_per_chain) || !is.numeric(n_iter_per_chain) || n_iter_per_chain <= 0) {
    stop("Error: 'n_iter_per_chain' parameter is missing or invalid. It must be a positive integer.")
  }
  if (n_chains > parallel::detectCores()) {
    stop("Error: 'n_chains' exceeds the number of available CPU cores.")
  }
  if (prev < 0 || prev > 1) {
    stop("Error: 'prev' must be between 0 and 1.")
  }

  # Check the length of var based on sex_specific
  if (sex_specific && length(var) != 8) {
    stop("Error: When 'sex_specific' is TRUE, 'var' must have exactly 8 elements.")
  } else if (!sex_specific && length(var) != 4) {
    stop("Error: When 'sex_specific' is FALSE, 'var' must have exactly 4 elements.")
  }

  # Check baseline_data structure
  if (sex_specific) {
    if (!is.data.frame(baseline_data) || nrow(baseline_data) != 94 || ncol(baseline_data) != 3) {
      stop("Error: 'baseline_data' must be a data frame with 94 rows and 3 columns (Age, Female, Male) when 'sex_specific' is TRUE.")
    }
  } else {
    if (!is.data.frame(baseline_data) && !is.vector(baseline_data)) {
      stop("Error: 'baseline_data' must be either a data frame with 94 rows and 1 column or a numeric vector when 'sex_specific' is FALSE.")
    }
    if (is.data.frame(baseline_data) && (nrow(baseline_data) != 94 || ncol(baseline_data) != 1)) {
      stop("Error: When 'baseline_data' is a data frame for non-sex-specific analysis, it must have 94 rows and 1 column.")
    }
    if (is.vector(baseline_data) && length(baseline_data) != 94) {
      stop("Error: When 'baseline_data' is a vector for non-sex-specific analysis, it must have exactly 94 elements.")
    }
  }

  # Create the seeds for the individual chains
  seeds <- sample.int(1000, n_chains)

  # Apply the transformation to adjust the format for the clipp package
  data <- do.call(rbind, lapply(pedigree, transformDF))

  # Create the prior distributions
  prop <- makePriors(
    data = distribution_data,
    sample_size = sample_size,
    ratio = ratio,
    prior_params = prior_params,
    risk_proportion = risk_proportion,
    baseline_data = baseline_data
  )

  cores <- parallel::detectCores()

  if (n_chains > cores) {
    stop("Error: 'n_chains exceeds the number of available CPU cores.")
  }
  cl <- parallel::makeCluster(n_chains)

  # Load required packages to the clusters
  parallel::clusterEvalQ(cl, {
    library(clipp)
    library(stats4)
    library(MASS)
    library(parallel)
    library(kinship2)
  })

  parallel::clusterExport(cl, c(
    "mhChain", "mhLogLikelihood_clipp", "mhLogLikelihood_clipp_noSex", "imputeUnaffectedAges",
    "calculate_weibull_parameters", "validate_weibull_parameters", "prior_params",
    "transformDF", "lik.fn", "lik_noSex", "mvrnorm", "var", "calculateEmpiricalDensity", "baseline_data",
    "seeds", "n_iter_per_chain", "burn_in", "imputeAges", "imputeAgesInit",
    "drawBaseline", "calculateNCPen", "drawEmpirical", "imp_interval",
    "data", "twins", "prop", "prev", "max_age", "BaselineNC", "median_max", "ncores",
    "remove_proband", "sex_specific"
  ), envir = environment())

  results <- parallel::parLapply(cl, 1:n_chains, function(i) {
    mhChain(
      seed = seeds[i],
      n_iter = n_iter_per_chain,
      burn_in = burn_in,
      chain_id = i,
      data = data,
      twins = twins,
      ncores = ncores,
      prior_distributions = prop,
      max_age = max_age,
      prev = prev,
      median_max = median_max,
      baseline_data = baseline_data,
      BaselineNC = BaselineNC,
      var = var,
      age_imputation = age_imputation,
      imp_interval = imp_interval,
      remove_proband = remove_proband,
      sex_specific = sex_specific
    )
  })

  # Check rejection rates and issue a warning if they are all above 90%
  all_high_rejections <- all(sapply(results, function(x) x$rejection_rate > 0.9))
  if (all_high_rejections) {
    warning("Low acceptance rate. Please consider running the chain longer.")
  }

  # Apply burn-in and thinning
  if (burn_in > 0) {
    results <- apply_burn_in(results, burn_in)
  }
  if (thinning_factor > 1) {
    results <- apply_thinning(results, thinning_factor)
  }

  # Select the appropriate combination chain function
  combine_function <- if (sex_specific) combine_chains else combine_chains_noSex

  # Select the appropriatesummary function
  summary_function <- if (sex_specific) generate_summary else generate_summary_noSex

  # Extract samples from the chains
  combined_chains <- combine_function(results)

  # Initialize variables
  output <- list()

  tryCatch(
    {
      if (rejection_rates) {
        # Generate rejection rates
        output$rejection_rates <- printRejectionRates(results)
      }

      if (summary_stats) {
        # Generate summary statistics
        output$summary_stats <- summary_function(combined_chains)
      }

      if (density_plots) {
        # Generate density plots
        output$density_plots <- generate_density_plots(combined_chains)
      }

      if (plot_trace) {
        # Generate trace plots
        output$trace_plots <- plot_trace(results, n_chains)
      }

      if (penetrance_plot) {
        # Generate penetrance plot
        output$penetrance_plot <- plot_penetrance(combined_chains, prob = probCI, max_age = max_age)
      }

      if (penetrance_plot_pdf) {
        # Generate PDF plots
        output$penetrance_plot_pdf <- plot_pdf(combined_chains, prob = probCI, max_age = max_age, sex = "NA")
      }

      if (plot_loglikelihood) {
        output$loglikelihood_plots <- plot_loglikelihood(results, n_chains)
      }

      if (plot_acf) {
        output$acf_plots <- plot_acf(results, n_chains)
      }
    },
    error = function(e) {
      # Handle errors here
      message("An error occurred in the output display: ", e$message)
    }
  )

  output$combined_chains <- combined_chains
  output$results <- results
  output$data <- data

  parallel::stopCluster(cl)

  return(output)
}