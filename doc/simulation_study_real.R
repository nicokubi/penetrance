## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(penetrance)
library(ggplot2)
library(scales)

## -----------------------------------------------------------------------------
dat <- test_fam2

## ----eval=FALSE---------------------------------------------------------------
#  
#  # Set the random seed
#  set.seed(2024)
#  
#  # Set the prior
#  prior_params <- list(
#      asymptote = list(g1 = 1, g2 = 1),
#      threshold = list(min = 5, max = 30),
#      median = list(m1 = 2, m2 = 2),
#      first_quartile = list(q1 = 6, q2 = 3)
#  )
#  
#  # Set the allele frequency for MLH1 based on PanelPRO Database
#  af_MLH1 <- 0.0004453125
#  
#  # We use the default baseline (non-carrier) penetrance
#  print(baseline_data_default)
#  
#  # We run the estimation procedure with one chain and 20k iterations
#  out_sim <- penetrance(
#      pedigree  = dat, twins = NULL, n_chains = 1, n_iter_per_chain = 20000,
#      ncores = 2, baseline_data = baseline_data_default , af = af_MLH1,
#      prior_params = prior_params, burn_in = 0.1, median_max = TRUE,
#      ageImputation = FALSE, removeProband = FALSE
#  )
#  

