################################################################################
#	Script-file:   contour_plot.R
# Aim:           Generate data for contour plot
# Project:       Capstone Project 2025
#
#	Data used:	   none
#	Data created:  data
# Student:       Zimo
################################################################################

# Clean the session
rm(list=ls()) 

# Change working dictionary
setwd("~/Documents/Biostatistics/project/R code")

library(tidyverse)
library(akima)
library(ggplot2)
library(patchwork)

# Simulate case counts function
simulate_case_counts <- function(R_t, W_s, I_0 = 10) {
  # R_t: A vector of reproduction numbers over time (length = number of days)
  # W_s: A vector representing the GI distribution
  # I_0: Initial number of cases (default is 10)
  
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
  
  # Return the simulated case counts
  return(I)
}

################################################################################
# Log-likelihood function
log_likelihood <- function(I, R_t, W_s) {
  # I: Daily observed case counts (length N)
  # R_t: Reproduction numbers (length N)
  # W_s: Generation interval probabilities (length <= N)
  
  n_days <- length(I) # Number of days
  E_I <- numeric(n_days)
  E_I[1] <- I[1] # Start with I[1] cases
  ll <- I[1] * log(E_I[1]) - E_I[1] # Initial value for log-likelihood
  
  # Calculate expected cases lambda_t
  for (t in 2:n_days) {
    
    E_I[t] <- 0
    for (s in 1:(t-1)) {
      
      # End if exceeds the length of W_s 
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

################################################################################
# Generate GI and Rt function
sim_GI_Rt <- function(GI_mean, GI_sd, Rt_mean, Rt_sd, Rt_length) {
  
  # Log-normal parameters for GI
  sdlog <- sqrt(log(1 + (GI_sd^2) / (GI_mean^2)))
  meanlog <- log(GI_mean) - (sdlog^2) / 2
  
  # Gamma parameter for reproduction number
  shape <- (Rt_mean^2) / (Rt_sd^2)
  scale <- (Rt_sd^2) / Rt_mean
  
  # Generate GI distribution
  W_s_raw <- dlnorm(1:50, meanlog, sdlog)
  W_s_raw <- W_s_raw[W_s_raw  > 1e-5]
  W_s <- W_s_raw/sum(W_s_raw) # standardise to make the sum equal to 1
  
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
# Generate log-likelihood data set function
ll_plot_data <- function(N = 100000, case_count_days = 50, I = case_count,
                         Rt_lower = 0.5, Rt_upper = 5, Rt_sd = 0.001){
  # N: the number of observations of mean, sd, and log-likelihood
  # case_count_days: the length of the simulated case count data
  # I: Simulated series of case count
 
   # generate a series of random number as mean and sd
  data <- tibble(
    GI_mean = runif(N, min = 3, max = 7), # GI mean range
    GI_sd = rep(1, N), # GI sd set to be 1
    Rt_mean = runif(N, min = Rt_lower, max = Rt_upper), # Rt mean range
    Rt_sd = rep(Rt_sd, N), # Rt sd set to be 0.001
    Rt_length = rep(case_count_days, N)
  )
  
  # generate the GI and Rt based on the simulated mean and sd
  results <- apply(data, 1, function(row) {
    do.call(sim_GI_Rt, as.list(row))
  })
  
  # combine the mean, sd, and GI, Rt into a data set 
  results_df <- map_dfr(results, as_tibble)
  final_data <- bind_cols(data, results_df)
  
  # Calculate the log-likelihood
  for (i in 1:N) {
    final_data$ll[i] <- log_likelihood(I, final_data$Rt[[i]], final_data$Ws[[i]])
  }
  
  return(final_data)
}

################################################################################
# Single plot
# Set parameters
parameters <- data.frame(
  GI_mean = 5,
  GI_sd = 1,
  Rt_mean = 2,
  Rt_sd = 0.001,
  Rt_lower = 0.5,
  Rt_upper = 3.5,
  days = 20
  )

# Simulate Rt and Ws
data <- sim_GI_Rt(GI_mean = parameters$GI_mean, GI_sd = parameters$GI_sd,  
                  Rt_mean = parameters$Rt_mean, Rt_sd = parameters$Rt_sd, 
                  Rt_length = parameters$days)

# Simulate case count data
case_count <- simulate_case_counts(data$Rt[[1]], data$Ws[[1]], I_0 = 10)

# Calculate log-likelihood
data <- ll_plot_data(N = 100000, case_count_days = parameters$days, I = case_count,
                     Rt_lower = parameters$Rt_lower, Rt_upper = parameters$Rt_upper, 
                     Rt_sd = parameters$Rt_sd)
data <- data[data$ll  >= -100,] # zoom in for log-likelihood < 100

# Generate regular grid data for contour plot
reg_data <- interp(
  x = data$GI_mean,   
  y = data$Rt_mean,   
  z = data$ll,  
  xo = seq(min(data$GI_mean), max(data$GI_mean), length = 100),  # X grid points
  yo = seq(min(data$Rt_mean), max(data$Rt_mean), length = 100),  # Y grid points
  linear = TRUE          # Use linear interpolation
)

# Convert to data frame for ggplot
reg_df <- data.frame(
  GI_mean = rep(reg_data$x, each = length(reg_data$y)),
  Rt_mean = rep(reg_data$y, length(reg_data$x)),
  log_like = as.vector(reg_data$z)
)

# Generate contour plot
ggplot(reg_df, aes(x = GI_mean, y = Rt_mean, z = log_like)) +
  geom_contour_filled() +  # Filled contours
  geom_contour() +  # Add contour lines
  scale_fill_viridis_d(name = "Log-Likelihood") +   # Color scale
  labs(title = "Log-Likelihood Contour Plot",
    x = "Generation Interval Mean (days)",
    y = "Reproduction Number (R)",
    caption = paste0("Simulated Case Count: ", toString(case_count))) +
  theme(plot.caption = element_text(hjust = 0)) +
  scale_x_continuous(limits = c(2, 7)) +
  geom_point(data=parameters, aes(x = GI_mean, y = Rt_mean), 
             size = 5, shape = 10, color = "#FF3300", inherit.aes = FALSE) +
  geom_text(data=parameters, aes(x = GI_mean, y = Rt_mean), 
            label=paste0("(",parameters$GI_mean,",",parameters$Rt_mean,")"), 
            nudge_y = 0.13, color = "#FF3300", inherit.aes = FALSE) +
  theme_minimal()

################################################################################
# Multiple plots
# Create an empty list for multiple plots
plot_list <- list()

# Set parameters
for (i in seq(1.5, 3, 0.5)){
  # Set parameters
  parameters <- data.frame(
    GI_mean = 5,
    GI_sd = 1,
    Rt_mean = i,
    Rt_sd = 0.001,
    Rt_lower = 1.5,
    Rt_upper = 3,
    days = 20
  )
  
  # Simulate Rt and Ws
  data <- sim_GI_Rt(GI_mean = parameters$GI_mean, GI_sd = parameters$GI_sd,  
                    Rt_mean = parameters$Rt_mean, Rt_sd = parameters$Rt_sd, 
                    Rt_length = parameters$days)
  
  # Simulate case count data
  case_count <- simulate_case_counts(data$Rt[[1]], data$Ws[[1]], I_0 = 10)
  
  # Calculate log-likelihood
  data <- ll_plot_data(N = 100000, case_count_days = parameters$days, I = case_count,
                       Rt_lower = parameters$Rt_lower, Rt_upper = parameters$Rt_upper, 
                       Rt_sd = parameters$Rt_sd)
  data <- data[data$ll  >= -100,] # zoom in for log-likelihood < 100
  
  # Generate regular grid data for contour plot
  reg_data <- interp(
    x = data$GI_mean,   
    y = data$Rt_mean,   
    z = data$ll,  
    xo = seq(min(data$GI_mean), max(data$GI_mean), length = 50),  # X grid points
    yo = seq(min(data$Rt_mean), max(data$Rt_mean), length = 50),  # Y grid points
    linear = TRUE          # Use linear interpolation
  )
  
  # Convert to data frame for ggplot
  reg_df <- data.frame(
    GI_mean = rep(reg_data$x, each = length(reg_data$y)),
    Rt_mean = rep(reg_data$y, length(reg_data$x)),
    log_like = as.vector(reg_data$z)
  )
  
  # Generate contour plot
  p <- ggplot(reg_df, aes(x = GI_mean, y = Rt_mean, z = log_like)) +
    geom_contour_filled() +  # Filled contours
    #geom_contour() +  # Add contour lines
    scale_fill_viridis_d(name = "Log-Likelihood") +  # Color scale
    labs(title = "Log-Likelihood Contour Plot",
         x = "Generation Interval Mean (days)",
         y = "Reproduction Number (R)",
         caption = paste0("Simulated Case Count: ", toString(case_count))) +
    theme(plot.caption = element_text(hjust = 0)) +
    scale_x_continuous(limits = c(2, 7)) +
    scale_y_continuous(limits = c(parameters$Rt_lower, parameters$Rt_upper)) +
    geom_point(data=parameters, aes(x = GI_mean, y = as.numeric(Rt_mean)), 
               size = 5, shape = 10, color = "#FF3300", inherit.aes = FALSE) +
    geom_text(data=parameters, aes(x = GI_mean, y = as.numeric(Rt_mean)), 
              label=paste0("(",parameters$GI_mean,",",parameters$Rt_mean,")"), 
              nudge_y = 0.13, color = "#FF3300", inherit.aes = FALSE) +
    theme_minimal()
  
  plot_list[[i*2]] <- p
}

# Combine the plots together
plot_list <- discard(plot_list, is.null)
combined_plot <- wrap_plots(plot_list, ncol = 2) # 2 columns
print(combined_plot)

# check GI shape and Rt length
w <- data.frame(w = data$Ws[[2908]], i = seq(1,length(data$Ws[[2908]]),1))
ggplot(w, aes(i, w)) + geom_line()
length(data$Rt[[2908]])
