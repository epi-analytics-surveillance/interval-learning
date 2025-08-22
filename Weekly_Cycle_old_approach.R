################################################################################
#	Script-file:   weekly_cycle_simulation.R
# Aim:           Function to simulate case counts with weekly pattern
# Project:       Capstone Project 2025
#
#	Data used:	   none
#	Data created:  none
# Student:       Zimo
################################################################################
# Clean the session
rm(list=ls()) 

# Change working dictionary
setwd("~/Documents/Biostatistics/project/R code")

library(dplyr)
library(ggplot2)
library(patchwork)

simulate_case_counts <- function(R_t, W_s, I_0 = 10, cycle_pattern = 0) {
  # R_t: A vector of reproduction numbers over time (length = number of days)
  # W_s: A vector representing the GI distribution
  # I_0: Initial number of cases (default is 10)
  # weekly_pattern: 0 - no pattern; other numbers - number of days in a cycle
  
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
  if (cycle_pattern != 0) {
    
    # 7 days as a week
    day_of_week <- (1:n_days - 1) %% cycle_pattern
    
    # Define pattern
    multiplier <- case_when(
      day_of_week %in% c(0, 1) ~ 0.75,  # Day 0 & 1: 0.75x (e.g., weekend)
      day_of_week == 2          ~ 1.5,  # Day 2:     1.3x  (e.g., Monday surge)
      TRUE                     ~ 1      # Other:     1x    (no adjustment)
    )
    
    # Multiply case counts
    I <- round(I * multiplier)
  }
  
  # Return the simulated case counts
  return(I)
}

################################################################################
# Generate GI and Rt function
sim_GI_Rt <- function(GI_mean, GI_sd, Rt_mean, Rt_sd, Rt_length) {
  
  # Log-normal parameters for GI
  sdlog <- sqrt(log(1 + (GI_sd^2) / (GI_mean^2)))
  meanlog <- log((GI_mean^2) / sqrt(GI_mean^2 + GI_sd^2))
  
  # Gamma parameter for reproduction number
  shape <- (Rt_mean^2) / (Rt_sd^2)
  scale <- (Rt_sd^2) / Rt_mean
  
  # Generate GI distribution
  W_s_raw <- dlnorm(1:50, meanlog, sdlog)
  W_s <- W_s_raw/sum(W_s_raw) # standardise to make the sum equal to 1
  W_s <- W_s[W_s  > 1e-20]
  
  # Generate simulated reproduction number
  R_t <- rgamma(Rt_length, shape = shape, scale = scale)
  
  # Return the result in list
  return(list(
    GI_sdlog = sdlog,
    GI_meanlog = meanlog,
    Rt_shape = shape,
    Rt_scale = scale,
    Ws = list(W_s),
    Rt = list(R_t)
  ))
}
################################################################################
# Log-likelihood function
log_likelihood <- function(I, R_t, W_s) {
  # I: Daily observed case counts (length N)
  # R_t: Reproduction numbers (length N)
  # W_s: Generation interval probabilities (length <= N)
  
  n_days <- length(I) # Number of days == length(R_t)
  E_I <- numeric(n_days)
  E_I[1] <- I[1] # Start with I[1] cases
  ll <- I[1] * log(E_I[1]) - E_I[1] - lgamma(I[1] + 1) # Initial value for log-likelihood
  
  # Calculate expected cases lambda_t
  for (t in 2:n_days) {
    
    E_I[t] <- 0
    for (s in 1:(t-1)) {
      
      # End if exceeds the length of W_s 
      if (s > length(W_s)) {
        break
      }
      
      E_I[t] <- E_I[t] + R_t[t - s] * W_s[s] * I[t - s]
    }
    
    
    # Calculate Poisson log-likelihood
    ll <- ll + I[t] * log(E_I[t]) - E_I[t] - lgamma(I[t] + 1)
  }
  
  return(ll)
}

################################################################################
# simulate case count
Rt <- rep(1.1, 200)
Ws <- sim_GI_Rt(GI_mean = 5, 
                GI_sd = 2, 
                0, 0, 0)$Ws[[1]]
case_count_whole <- simulate_case_counts(Rt, Ws, I_0 = 1000, cycle_pattern = 7)
plot(case_count_whole, type = "l")

# group the case count into 7 days
n <- length(case_count_whole)  
num_repeats <- 7
total_weeks <- ceiling(n / num_repeats)
sequence <- rep(1:total_weeks, each = num_repeats)[1:n]
series_data <- data.frame(weeks = sequence, case = case_count_whole)
groups <- split(series_data$case, as.factor(series_data$weeks))

Rt_sequence <- data.frame()
contour_plot <- list()

# Infer Rt and GI
for (i in 1:total_weeks) {
  
  # case count
  case_count7 = groups[[i]]
  
  # Set GI parameters and Rt mean
  parameters_data <- expand.grid(GI_mean = seq(1, 10, 1), 
                                 GI_sd = seq(0.5, 5, 0.5), 
                                 Rt_mean = seq(0.5, 3, 0.2))
  
  # Simulate GI and Rt, and compute log-likelihood for every 7 days
  for (j in 1:nrow(parameters_data)) {
    
    Result <- sim_GI_Rt(parameters_data$GI_mean[j], parameters_data$GI_sd[j], 
                        parameters_data$Rt_mean[j], 0.01, 
                        length(groups[[i]]))
    parameters_data$GI[j] <- Result$Ws
    parameters_data$Rt[j] <- Result$Rt
    parameters_data$ll[j] <- log_likelihood(case_count7, 
                                            Result$Rt[[1]], 
                                            Result$Ws[[1]])
  }
  
  # Select the best fit GI sd for each pair of GI and Rt mean
  optimal_parameters <- parameters_data %>%
    group_by(GI_mean, Rt_mean) %>%
    slice(which.max(ll)) %>%
    ungroup() %>%
    dplyr::select(GI_mean, GI_sd, GI, Rt_mean, Rt, ll) %>%
    arrange(GI_mean)
  
  # Generate contour plot
  p <- ggplot(optimal_parameters, aes(x = GI_mean, y = Rt_mean, z = ll)) +
    geom_contour_filled(bins=20) +  # Filled contours
    geom_contour() +  # Add contour lines
    scale_fill_viridis_d(name = "Log-Likelihood") +   # Color scale
    labs(title = paste0("Log-Likelihood Contour Plot: Week ",i),
         x = "Generation Interval Mean (days)",
         y = "Reproduction Number") + 
    theme_minimal()
  
  contour_plot[[i]] <- p
  
  # Extract the row with maximum log-likelihood
  best_fit <- parameters_data[parameters_data$ll == max(parameters_data$ll),]
  
  # Append to Rt_data (ensure column names match)
  Rt_sequence <- rbind(Rt_sequence, best_fit)
}

# View contour plot for the first week
contour_plot[[1]] / contour_plot[[19]] |
  contour_plot[[9]] / contour_plot[[29]]
