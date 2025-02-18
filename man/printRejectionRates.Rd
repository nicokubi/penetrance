\name{printRejectionRates}
\alias{printRejectionRates}
\title{Print MCMC Rejection Rates}
\usage{
printRejectionRates(results)
}
\arguments{
  \item{results}{A list of MCMC chain results.}
}
\value{
  A numeric vector containing the rejection rates for each chain.
}
\description{
  Calculates and prints the rejection rates for each MCMC chain.
}
\examples{
  # Create example chain results
  results <- list(
    list(rejection_rate = 0.3),
    list(rejection_rate = 0.4)
  )
  
  # Print rejection rates
  rates <- printRejectionRates(results)
}