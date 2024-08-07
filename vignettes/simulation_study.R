## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(penetrance)
library(ggplot2)

## -----------------------------------------------------------------------------
# Create generating_penetrance data frame
age <- 1:94

# Calculate Weibull distribution for Females
alpha <- 2 
beta  <- 50 
gamma <- 0.6
delta <- 15
penetrance.mod.f <- dweibull(age - delta, alpha, beta) * gamma

# Calculate Weibull distribution for Males
alpha <- 2 
beta  <- 50 
gamma <- 0.6
delta <- 30
penetrance.mod.m <- dweibull(age - delta, alpha, beta) * gamma

generating_penetrance <- data.frame(
    Age = age,
    Female = penetrance.mod.f,
    Male = penetrance.mod.m
)

## -----------------------------------------------------------------------------

dat <- test_fam


