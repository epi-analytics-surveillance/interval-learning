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
# Generate GI function
sim_GI <- function(GI_mean, GI_sd) {
  
  # Log-normal parameters for GI
  sdlog <- sqrt(log(1 + (GI_sd^2) / (GI_mean^2)))
  meanlog <- log((GI_mean^2) / sqrt(GI_mean^2 + GI_sd^2))
  
  # Generate GI distribution
  W_s <- dlnorm(1:50, meanlog, sdlog)
  W_s <- W_s/sum(W_s) # standardise to make the sum equal to 1
  
  # Return the result in list
  return(list(
    GI_sdlog = sdlog,
    GI_meanlog = meanlog,
    Ws = list(W_s)
  ))
}
################################################################################
# Log-likelihood function
log_likelihood_daily <- function(t, I, Rt, Ws) {
  # t: day index
  # I: Daily observed case counts 
  # R_t: Reproduction numbers
  # W_s: Generation interval probabilities 
  
  # Set the infectiousness or namely the potential
  lambda_t <- 0
  
  # Calculate lambda_t
  for (s in 1:min(length(Ws), t-1)) {
    if ((t - s) >= 1) {
      lambda_t <- lambda_t + Ws[s] * I[t - s]
    }
  }
  
  # Calculate expected cases on day t
  if (t == 1) {
    EI <- I[1]
  } else {
    EI <- Rt * lambda_t
  }
  
  # Handle extremely small EI
  if (EI < 1e-100) {
    if (I[t] == 0) {
      return(-EI)  # log(P(0)) ≈ -EI when EI is small
    } else {
      return(-1e100)  # Virtually impossible
    }
  }
  
  # Calculate log-likelihood on day t
  ll_t <- I[t] * log(EI) - EI - lgamma(I[t] + 1)
  
  return(ll_t)
}

################################################################################
# simulate case count
sim_p <- data.frame(
  GI_mean = 6,
  GI_sd = 1.5,
  Rt_mean = 1.2,
  cycle_pattern = 4
)

Rt <- rep(sim_p$Rt_mean, 200)
Ws <- sim_GI(sim_p$GI_mean, sim_p$GI_sd)$Ws[[1]]
case_count_whole <- simulate_case_counts(Rt, Ws, I_0 = 1000, sim_p$cycle_pattern)
plot(case_count_whole, type = "l")

# group the case count into 7 days
n <- length(case_count_whole)  
num_repeats <- 7
total_weeks <- ceiling(n / num_repeats)
sequence <- rep(1:total_weeks, each = num_repeats)[1:n]
series_data <- data.frame(weeks = sequence, case = case_count_whole)
groups <- split(series_data$case, as.factor(series_data$weeks))

# Pre-set to store information from the for loop
Rt_sequence <- data.frame()
contour_plot <- list()

# Infer Rt and GI
for (i in 1:total_weeks) {
  
  # Set GI parameters and Rt mean
  parameters_data <- expand.grid(GI_mean = seq(1, 10, 1), 
                                 GI_sd = seq(0.5, 5, 0.5), 
                                 Rt_mean = seq(0.5, 3, 0.2))
  
  # Simulate GI and compute log-likelihood for every 7 days
  for (j in 1:nrow(parameters_data)) {
    
    # Simulate GI distribution for each pair of GI mean and sd
    Result <- sim_GI(parameters_data$GI_mean[j], parameters_data$GI_sd[j])
    parameters_data$GI[j] <- Result$Ws
    
    daily_ll <- numeric()
    
    # Calculate daily log-likelihood
    for (k in 1:length(groups[[i]])) {
      t <- 7*(i-1) + k  # start from the first day of each group
      daily_ll[k] <- log_likelihood_daily(t, case_count_whole, 
                                          parameters_data$Rt_mean[j], Result$Ws[[1]])
    }
    
    # Sum up daily log-likelihood for the week
    parameters_data$ll[j] <- sum(daily_ll)
  }
  
  # Select the best fit GI sd for each pair of GI and Rt mean
  optimal_parameters <- parameters_data %>%
    group_by(GI_mean, Rt_mean) %>%
    slice(which.max(ll)) %>%
    ungroup() %>%
    dplyr::select(Rt_mean, GI_mean, GI_sd, GI, ll) %>%
    arrange(GI_mean)
  
  # Zoom in the plot
  threshold <- quantile(optimal_parameters$ll, probs = 0.5)
  optimal_parameters <- optimal_parameters[optimal_parameters$ll >= threshold, ]
  
  # Generate contour plot
  p <- ggplot(optimal_parameters, aes(x = GI_mean, y = Rt_mean, z = ll)) +
    geom_contour_filled(bins=20) +  # Filled contours
    geom_contour() +  # Add contour lines
    scale_fill_viridis_d(name = "Log-Likelihood") +   # Color scale
    labs(title = paste0("Log-Likelihood Contour Plot: Week ",i, 
                        "; Cycle: ", sim_p$cycle_pattern, " days"),
         x = "Generation Interval Mean (days)",
         y = "Reproduction Number") + 
    geom_point(data=sim_p, aes(x = GI_mean, y = as.numeric(Rt_mean)), 
               size = 5, shape = 10, color = "#FF3300", inherit.aes = FALSE) +
    geom_text(data=sim_p, aes(x = GI_mean, y = as.numeric(Rt_mean)), 
              label=paste0("(",sim_p$GI_mean,",",sim_p$Rt_mean,")"), 
              nudge_y = 0.13, color = "#FF3300", inherit.aes = FALSE) +
    theme_minimal()
  # store the plot for each week
  contour_plot[[i]] <- p
  
  # Extract the row with maximum log-likelihood
  best_fit <- parameters_data[parameters_data$ll == max(parameters_data$ll),]
  # check error
  if (nrow(best_fit) > 1) {
    print(best_fit)
    break
  }
  
  # Append to Rt_data (ensure column names match)
  Rt_sequence <- rbind(Rt_sequence, best_fit)
}

# View contour plot for the first week
contour_plot[[1]] / contour_plot[[19]] |
  contour_plot[[9]] / contour_plot[[29]]

# Check the best fit parameters for each week
View(Rt_sequence)
