####################################################################################
#	Script-file:   log-likelihood.R
# Aim:           Simulate incidence data from Poisson distribution
# Project:       Capstone Project 2025
#
#	Data used:	   none
#	Data created:  none
# Student:       Zimo
####################################################################################

# Clean the session
rm(list=ls()) 

# Change working dictionary
setwd("~/Documents/Biostatistics/project/R code")

# Calculate log-likelihood
log_likelihood <- function(I, R_t, W_s) {
  # I: Daily observed case counts (length N)
  # R_t: Reproduction numbers (length N)
  # w: Generation interval probabilities (length <= N)
  
  n_days <- length(I) # Number of days
  E_I <- numeric(n_days)
  E_I[1] <- I[1] # Start with I[1] cases
  ll <- I[1] * log(E_I[1]) - E_I[1] # Initial value for log-likelihood
  
  # Calculate expected cases lambda_t
  for (t in 2:n_days) {
    
    E_I[t] <- 0
    for (s in 1:(t-1)) {
      
      # End if the exceeds length of W_s 
      if (s > length(W_s)) {
        break
      }
      
      E_I[t] <- E_I[t] + R_t[s] * I[t - s] * W_s[s]
    }
    
    # Calculate Poisson log-likelihood
    ll <- ll + I[t] * log(E_I[t]) - E_I[t] - lgamma(I[t] + 1)
  }
  
  return(ll)
}
