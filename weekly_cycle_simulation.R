################################################################################
#	Script-file:   weekly_cycle_simulation.R
# Aim:           Function to simulate case counts with weekly pattern
# Project:       Capstone Project 2025
#
#	Data used:	   none
#	Data created:  none
# Student:       Zimo
################################################################################
simulate_case_counts <- function(R_t, W_s, I_0 = 10, weekly_pattern = TRUE) {
  # R_t: A vector of reproduction numbers over time (length = number of days)
  # W_s: A vector representing the GI distribution
  # I_0: Initial number of cases (default is 10)
  # weekly_pattern: Whether to apply weekly pattern (default is TRUE)
  
  # Number of days
  n_days <- length(R_t)
  
  # Initialize the case count vector
  I <- numeric(n_days)
  I[1] <- I_0  # Start with the initial number of cases on day 1
  
  # Simulate cases for subsequent days
  for (t in 2:n_days) {
    # Calculate the total infection potential Lambda_t
    lambda_t <- 0
    for (s in 1:(t-1)) {
      
      # End if exceeds length of W_s
      if (s > length(W_s)) {
        break
      }
      
      # Calculate the expected number of cases on day t
      lambda_t <- lambda_t + I[t - s] * W_s[s] * R_t[t - s]
    }
    
    # Simulate new cases from a Poisson distribution
    I[t] <- rpois(1, lambda = lambda_t)
  }
  
  # Apply weekly pattern if requested
  if (weekly_pattern) {
    
    # 7 days as a week
    day_of_week <- (1:n_days - 1) %% 7  # 0,1,2,...,6 (repeat every 7 days)
    
    # Define pattern
      multiplier <- case_when(
        day_of_week %in% c(5, 6) ~ 0.75,  # Days 0,1: 0.75 (e.g., weekend)
        day_of_week == 1          ~ 1.3,  # Day 2: 1.3x    (e.g., Monday surge)
        TRUE                     ~ 1      # Days 3-6: 1x   (no adjustment)
      )
    
    # Multiply case counts
    I <- round(I * multiplier)
  }
  
  # Return the simulated case counts
  return(I)
}

Rt <- rep(1.3, 40)
Ws <- c(0.1, 0.2, 0.4, 0.2, 0.1)
case <- simulate_case_counts(Rt, Ws, 10, TRUE)
plot(seq(1, 40, 1), case)
