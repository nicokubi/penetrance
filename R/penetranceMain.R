#' Bayesian Inference using Independent Metropolis-Hastings for Penetrance Estimation
#'
#' This function implements the Independent Metropolis-Hastings algorithm for Bayesian
#' penetrance estimation of cancer risk. It utilizes parallel computing to run multiple
#' chains and provides various options for analyzing and visualizing the results.
#'
#' @param pedigree A list of data frames, where each data frame represents a single pedigree and contains the following columns:
#'   - `PedigreeID`: A numeric or character identifier for the family/pedigree. Must be consistent for all members of the same family within a data frame.
#'   - `ID`: A unique numeric or character identifier for each individual within their respective pedigree data frame.
#'   - `Sex`: An integer representing biological sex: `0` for female, `1` for male. Use `NA` for unknown sex.
#'   - `MotherID`: The `ID` of the individual's mother. Should correspond to an `ID` within the same pedigree data frame or be `NA` if the mother is not in the pedigree (founder).
#'   - `FatherID`: The `ID` of the individual's father. Should correspond to an `ID` within the same pedigree data frame or be `NA` if the father is not in the pedigree (founder).
#'   - `isProband`: An integer indicating if the individual is a proband: `1` for proband, `0` otherwise.
#'   - `CurAge`: An integer representing the age of censoring. This is the current age if the individual is alive, or the age at death if deceased. Must be between `1` and `max_age`. Use `NA` for unknown ages (but note this may affect analysis or require imputation).
#'   - `isAff`: An integer indicating the affection status for the cancer of interest: `1` if diagnosed, `0` if unaffected. Use `NA` for unknown status.
#'   - `Age`: An integer representing the age at cancer diagnosis. Should be `NA` if `isAff` is `0` or `NA`. Must be between `1` and `max_age`, and less than or equal to `CurAge`. Use `NA` for unknown diagnosis age (but note this may affect analysis or require imputation).
#'   - `Geno`: An integer representing the germline genetic test result: `1` for carrier (positive), `0` for non-carrier (negative). Use `NA` for unknown or untested individuals.
#' @param twins A list specifying identical twins or triplets in the family. Each element of the list should be a vector containing the `ID`s of the identical siblings within a pedigree. For example: `list(c("ID1", "ID2"), c("ID3", "ID4", "ID5"))`. Default is `NULL`.
#' @param n_chains Integer, the number of chains for parallel computation. Default is 1.
#' @param n_iter_per_chain Integer, the number of iterations for each chain. Default is 10000.
#' @param ncores Integer, the number of cores for parallel computation. Default is 6.
#' @param baseline_data Data providing the absolute age-specific baseline risk (probability) of developing the cancer in the general population (e.g., from SEER database).
#'                      All probability values must be between 0 and 1.
#'                      - If `sex_specific = TRUE` (default): A data frame with columns 'Male' and 'Female', where each column contains the age-specific probabilities for that sex. The number of rows should ideally correspond to `max_age`.
#'                      - If `sex_specific = FALSE`: A numeric vector or a single-column data frame containing the age-specific probabilities for the combined population. The length (or number of rows) should ideally correspond to `max_age`.
#'                      Default data is provided for Colorectal cancer from SEER (up to age 94). If the number of rows/length does not match `max_age`, the data will be truncated or extended with the last value.
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
#'   Geno = c(NA, NA, 1, NA)
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
  if (!all(sapply(pedigree, is.data.frame))) {
      stop("Error: Each element in the 'pedigree' list must be a data frame.")
  }

  # Validate max_age
  if (!is.numeric(max_age) || length(max_age) != 1 || max_age <= 0 || max_age > 150 || floor(max_age) != max_age) {
    stop("Error: 'max_age' must be a single positive integer not exceeding 150.")
  }

  # Validate logical flags
  logical_params <- list(remove_proband = remove_proband, age_imputation = age_imputation, 
                         median_max = median_max, BaselineNC = BaselineNC, 
                         summary_stats = summary_stats, rejection_rates = rejection_rates, 
                         density_plots = density_plots, plot_trace = plot_trace, 
                         penetrance_plot = penetrance_plot, penetrance_plot_pdf = penetrance_plot_pdf,
                         plot_loglikelihood = plot_loglikelihood, plot_acf = plot_acf,
                         sex_specific = sex_specific)
  for (param_name in names(logical_params)) {
    param_value <- logical_params[[param_name]]
    if (!is.logical(param_value) || length(param_value) != 1 || is.na(param_value)) {
      stop(paste("Error: '", param_name, "' must be a single TRUE or FALSE value.", sep=""))
    }
  }

  # Validate other numeric/integer parameters
  if (!is.numeric(n_chains) || length(n_chains) != 1 || n_chains <= 0 || floor(n_chains) != n_chains) {
    stop("Error: 'n_chains' parameter must be a single positive integer.")
  }
  if (!is.numeric(n_iter_per_chain) || length(n_iter_per_chain) != 1 || n_iter_per_chain <= 0 || floor(n_iter_per_chain) != n_iter_per_chain) {
    stop("Error: 'n_iter_per_chain' parameter must be a single positive integer.")
  }
  if (!is.numeric(ncores) || length(ncores) != 1 || ncores <= 0 || floor(ncores) != ncores) {
    stop("Error: 'ncores' parameter must be a single positive integer.")
  }
   detected_cores <- parallel::detectCores()
  if (n_chains > detected_cores) {
    warning(paste("'n_chains' (", n_chains, ") exceeds the number of available CPU cores (", detected_cores, ").", sep=""))
    # Consider stopping if n_chains > detected_cores, but warning allows override.
    # stop("Error: 'n_chains' exceeds the number of available CPU cores.") 
  }
   if (ncores > detected_cores) {
    warning(paste("'ncores' (", ncores, ") exceeds the number of available CPU cores (", detected_cores, "). Using ", detected_cores, " cores instead.", sep=""))
    ncores <- detected_cores
  }
  if (!is.numeric(prev) || length(prev) != 1 || prev < 0 || prev > 1) {
    stop("Error: 'prev' must be a single numeric value between 0 and 1.")
  }
  if (!is.numeric(burn_in) || length(burn_in) != 1 || burn_in < 0 || burn_in >= 1) {
    stop("Error: 'burn_in' must be a single numeric value between 0 (inclusive) and 1 (exclusive).")
  }
  if (!is.numeric(thinning_factor) || length(thinning_factor) != 1 || thinning_factor <= 0 || floor(thinning_factor) != thinning_factor) {
    stop("Error: 'thinning_factor' must be a single positive integer.")
  }
   if (age_imputation && (!is.numeric(imp_interval) || length(imp_interval) != 1 || imp_interval <= 0 || floor(imp_interval) != imp_interval)) {
    stop("Error: 'imp_interval' must be a single positive integer when 'age_imputation' is TRUE.")
  }
  if (!is.numeric(probCI) || length(probCI) != 1 || probCI <= 0 || probCI >= 1) {
    stop("Error: 'probCI' must be a single numeric value between 0 (exclusive) and 1 (exclusive).")
  }

  # Validate var length based on sex_specific
  expected_var_length <- if (sex_specific) 8 else 4
  if (!is.numeric(var) || length(var) != expected_var_length) {
     stop(paste("Error: 'var' must be a numeric vector of length", expected_var_length, "when 'sex_specific' is", sex_specific))
  }
  if (any(var <= 0)) {
     stop("Error: All values in 'var' (proposal variances) must be positive.")
  }
  
  # Validate prior_params structure (basic check)
  if (!missing(prior_params) && !is.list(prior_params)) {
    stop("Error: 'prior_params' must be a list.")
  }
  # Add more specific checks for prior_params contents if needed, e.g. prior_params$shape, prior_params$scale

  # Validate pedigree data structure and content for each pedigree in the list
  required_columns <- c("PedigreeID", "ID", "Sex", "MotherID", "FatherID", "isProband", "CurAge", "isAff", "Age", "Geno")
  for (i in seq_along(pedigree)) {
    ped_df <- pedigree[[i]]
    ped_id_for_error <- unique(ped_df$PedigreeID)[1] # Get a representative PedigreeID for error messages

    if (!all(required_columns %in% colnames(ped_df))) {
      missing_cols <- setdiff(required_columns, colnames(ped_df))
      stop(paste("Error: Pedigree", ped_id_for_error, "(index", i, ") is missing required column(s):", paste(missing_cols, collapse=", ")))
    }

    # Check for NA/invalid values in critical columns
    critical_id_columns <- c("PedigreeID", "ID")
    for (col in critical_id_columns) {
      if (any(is.na(ped_df[[col]]))) {
        stop(paste("Error: NA values found in the '", col, "' column of pedigree ", ped_id_for_error, " (index ", i, ").", sep=""))
      }
    }
     if (any(duplicated(ped_df$ID))) {
         stop(paste("Error: Duplicated IDs found in pedigree ", ped_id_for_error, " (index ", i, "). IDs must be unique within a pedigree.", sep=""))
     }

    # Check Sex values
    if (!all(ped_df$Sex %in% c(0, 1, NA))) {
      stop(paste("Error: 'Sex' column in pedigree", ped_id_for_error, "(index", i, ") must only contain 0 (female), 1 (male), or NA (unknown)."))
    }

    # Check isProband values
    if (!all(ped_df$isProband %in% c(0, 1, NA))) { # Allow NA for flexibility, though 0/1 is standard
        stop(paste("Error: 'isProband' column in pedigree", ped_id_for_error, "(index", i, ") should only contain 0, 1, or NA."))
    }
     if (sum(ped_df$isProband == 1, na.rm = TRUE) == 0) {
        warning(paste("Warning: No proband (isProband=1) found in pedigree", ped_id_for_error, "(index", i, ")."))
    }


    # Check isAff values
    if (!all(ped_df$isAff %in% c(0, 1, NA))) {
      stop(paste("Error: 'isAff' column in pedigree", ped_id_for_error, "(index", i, ") must only contain 0 (unaffected), 1 (affected), or NA (unknown)."))
    }

    # Check geno values
    if (!all(ped_df$Geno %in% c(0, 1, NA))) {
      stop(paste("Error: 'Geno' column in pedigree", ped_id_for_error, "(index", i, ") must only contain 0 (negative/non-carrier), 1 (positive/carrier), or NA (unknown)."))
    }
    
    # Validate and attempt to coerce CurAge column to numeric
    original_na_curage <- sum(is.na(ped_df$CurAge))
    ped_df$CurAge <- suppressWarnings(as.numeric(as.character(ped_df$CurAge)))
    # Round to nearest integer
    ped_df$CurAge <- round(ped_df$CurAge)
    new_na_curage <- sum(is.na(ped_df$CurAge))
    if (new_na_curage > original_na_curage) {
        warning(paste("Warning: Non-numeric values found in 'CurAge' for pedigree", ped_id_for_error, "(index", i, ") were coerced to NA."))
    }
    # Check CurAge values are integers within the valid range [1, max_age]
    # Note: Rounding handles the integer requirement, just check the range.
    invalid_curage <- !is.na(ped_df$CurAge) & (ped_df$CurAge < 1 | ped_df$CurAge > max_age)
    if (any(invalid_curage)) {
        stop(paste("Error: 'CurAge' in pedigree", ped_id_for_error, "(index", i, ") contains values outside the valid range [1, ", max_age, "] after rounding. Check individuals: ", paste(ped_df$ID[invalid_curage], collapse=", ")))
    }
    
    # Validate and attempt to coerce Age column to numeric
    original_na_age <- sum(is.na(ped_df$Age))
    ped_df$Age <- suppressWarnings(as.numeric(as.character(ped_df$Age)))
    # Round to nearest integer
    ped_df$Age <- round(ped_df$Age)
    new_na_age <- sum(is.na(ped_df$Age))
     if (new_na_age > original_na_age) {
        warning(paste("Warning: Non-numeric values found in 'Age' (diagnosis age) for pedigree", ped_id_for_error, "(index", i, ") were coerced to NA."))
    }
    # Check Age values are integers within the valid range [1, max_age]
    # Note: Rounding handles the integer requirement, just check the range.
    invalid_age <- !is.na(ped_df$Age) & (ped_df$Age < 1 | ped_df$Age > max_age)
     if (any(invalid_age)) {
        stop(paste("Error: 'Age' (age of diagnosis) in pedigree", ped_id_for_error, "(index", i, ") contains values outside the valid range [1, ", max_age, "] after rounding. Check individuals: ", paste(ped_df$ID[invalid_age], collapse=", ")))
    }

    # Check consistency between Age, CurAge, and isAff (using rounded numeric values)
    # Note: Check if rounding causes Age > CurAge issues, although unlikely if original data was consistent.
    age_inconsistency <- !is.na(ped_df$Age) & !is.na(ped_df$CurAge) & (ped_df$Age > ped_df$CurAge)
     if (any(age_inconsistency)) {
        # It might be worth warning instead of stopping if rounding caused this.
        warning(paste("Warning: Age of diagnosis ('Age') may be greater than current age ('CurAge') after rounding for some individuals in pedigree", ped_id_for_error, "(index", i, "). Check individuals: ", paste(ped_df$ID[age_inconsistency], collapse=", ")))
        # Alternatively, stop as before:
        # stop(paste("Error: Age of diagnosis ('Age') is greater than current age ('CurAge') for some individuals in pedigree", ped_id_for_error, "(index", i, "). Check individuals: ", paste(ped_df$ID[age_inconsistency], collapse=", ")))
    }
    
     # Check MotherID and FatherID correspond to IDs within the same pedigree (optional but good practice)
    all_ids <- ped_df$ID
    invalid_mother_ids <- !is.na(ped_df$MotherID) & !(ped_df$MotherID %in% all_ids)
    invalid_father_ids <- !is.na(ped_df$FatherID) & !(ped_df$FatherID %in% all_ids)
    if (any(invalid_mother_ids)) {
       warning(paste("Warning: Some MotherIDs in pedigree", ped_id_for_error, "(index", i, ") do not correspond to existing IDs in the same pedigree. Ensure founders have NA parent IDs. Invalid MotherIDs for individuals:", paste(ped_df$ID[invalid_mother_ids], collapse=", ")))
    }
     if (any(invalid_father_ids)) {
       warning(paste("Warning: Some FatherIDs in pedigree", ped_id_for_error, "(index", i, ") do not correspond to existing IDs in the same pedigree. Ensure founders have NA parent IDs. Invalid FatherIDs for individuals:", paste(ped_df$ID[invalid_father_ids], collapse=", ")))
    }
    
    # Assign potentially modified ped_df back to the list to reflect numeric coercion
    pedigree[[i]] <- ped_df 
  }

  # Validate baseline_data structure and values
  if (sex_specific) {
    if (!is.data.frame(baseline_data)) {
      stop("Error: 'baseline_data' must be a data frame when 'sex_specific' is TRUE.")
    }
    required_bl_columns <- c("Male", "Female")
    if (!all(required_bl_columns %in% colnames(baseline_data))) {
      stop("Error: 'baseline_data' must have columns named 'Male' and 'Female' when 'sex_specific' is TRUE.")
    }
     if (!all(sapply(baseline_data[, required_bl_columns], is.numeric))) {
         stop("Error: 'Male' and 'Female' columns in 'baseline_data' must be numeric.")
     }
     if (any(baseline_data$Male < 0, na.rm = TRUE) || any(baseline_data$Male > 1, na.rm = TRUE) ||
         any(baseline_data$Female < 0, na.rm = TRUE) || any(baseline_data$Female > 1, na.rm = TRUE)) {
         stop("Error: Baseline probabilities in 'baseline_data' must be between 0 and 1.")
     }
    # Check if baseline_data matches max_age (warning handled below)
    data_rows <- nrow(baseline_data)
  } else { # Not sex_specific
    if (is.data.frame(baseline_data)) {
      if (ncol(baseline_data) != 1) {
        stop("Error: When 'baseline_data' is a data frame and 'sex_specific' is FALSE, it must have exactly one column.")
      }
      bl_vector <- baseline_data[[1]]
    } else if (is.vector(baseline_data) && is.numeric(baseline_data)) {
      bl_vector <- baseline_data
    } else {
      stop("Error: 'baseline_data' must be a numeric vector or a single-column data frame when 'sex_specific' is FALSE.")
    }
     if (any(bl_vector < 0, na.rm = TRUE) || any(bl_vector > 1, na.rm = TRUE)) {
         stop("Error: Baseline probabilities in 'baseline_data' must be between 0 and 1.")
     }
    data_rows <- length(bl_vector)
  }
  
   # Common baseline_data row check and warning
   if (data_rows < max_age) {
      warning(paste("Baseline data has fewer entries (", data_rows, ") than max_age (", max_age, "). Data will be extended by repeating the last value.", sep=""))
    } else if (data_rows > max_age) {
      warning(paste("Baseline data has more entries (", data_rows, ") than max_age (", max_age, "). Data will be truncated to the first ", max_age, " entries.", sep=""))
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
    stop("Error: 'n_chains' exceeds the number of available CPU cores.")
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