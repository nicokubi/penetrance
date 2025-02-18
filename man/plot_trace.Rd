\name{plot_trace}
\alias{plot_trace}
\title{Plot MCMC Trace Plots}
\usage{
plot_trace(results, n_chains)
}
\arguments{
  \item{results}{A list of MCMC chain results.}
  \item{n_chains}{The number of chains.}
}
\value{
  No return value, called for side effects. Creates trace plots for each parameter in each chain.
}
\description{
  Generates trace plots for visualizing MCMC chain convergence.
}
\examples{
  # Create example chain results
  results <- list(
    list(
      median_samples = rnorm(100),
      threshold_samples = runif(100),
      first_quartile_samples = rgamma(100, 2, 2),
      asymptote_samples = rbeta(100, 2, 2)
    )
  )
  
  # Plot traces for a single chain
  plot_trace(results, n_chains = 1)
}