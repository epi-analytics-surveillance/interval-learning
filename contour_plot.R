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
# transform into log-normal and gamma distribution parameters
log_normal_parameter <- function(GI_mean, GI_sd, Rt_mean , Rt_sd) {
  # Log-normal parameters for GI
  sdlog <- sqrt(log(1 + (GI_sd^2) / (GI_mean^2)))
  meanlog <- log(GI_mean) - (sdlog^2) / 2
  
  # Gamma parameter for reproduction number
  shape <- (Rt_mean^2) / (Rt_sd^2)
  scale <- (Rt_sd^2) / Rt_mean
  
  # Return the result in list
  return(list(
    GI_sdlog = sdlog,
    GI_meanlog = meanlog,
    Rt_shape = shape,
    Rt_scale = scale
  ))
}

################################################################################
# Generate GI and Rt
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

# data for 20 days
data <- sim_GI_Rt(GI_mean = 4.5, GI_sd = 1,  
                  Rt_mean = 2.75, Rt_sd = 0.001, Rt_length = 20)

# Simulate case count data
case_count <- simulate_case_counts(data$Rt[[1]], data$Ws[[1]], I_0 = 10)
################################################################################
# Generate log-likelihood data set
ll_plot_data <- function(N = 100000, case_count_days = 20, I = case_count){
  # N: the number of observations of mean, sd, and log-likelihood
  # case_count_days: the length of the simulated case count data
  # I: Simulated series of case count
 
   # generate a series of random number as mean and sd
  data <- tibble(
    GI_mean = runif(N, min = 2, max = 7), # GI mean range from 2-7 days
    GI_sd = runif(N, min = 1, max = 3), # GI sd range from 1-3 days
    Rt_mean = runif(N, min = 2, max = 7), # Rt mean range from 0.5-5
    Rt_sd = rep(0.001, N), # Rt sd set to be 0.001
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
    final_data$log_like[i] <- log_likelihood(I, final_data$Rt[[i]], final_data$Ws[[i]])
  }
  
  return(final_data)
}

data <- ll_plot_data()
View(data) # check the parameters, Rt, Ws, and log-likelihood

# Generate regular grid data for contour plot
reg_data <- interp(
  x = data$GI_mean,   
  y = data$Rt_mean,   
  z = data$log_like,  
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
ggplot(reg_df, aes(x = GI_mean, y = Rt_mean, z = log_like)) +
  geom_contour_filled() +  # Filled contours
  geom_contour() +  # Add contour lines
  scale_fill_viridis_d(name = "Log-Likelihood") +   # Color scale
  labs(
    title = "Log-Likelihood Contour Plot",
    x = "Generation Interval Mean (days)",
    y = "Reproduction Number (R)"
  ) +
  theme_minimal()


# Alternative: use filled.contour
filled.contour(x = reg_data$x,
               y = reg_data$y,
               z = reg_data$z,
               color.palette =
                 colorRampPalette(c("white", "blue")),
               xlab = "Generation Interval Mean (days)",
               ylab = "Reproduction Number (R)",
               main = "Log-Likelihood Contour Plot",
               key.title = title(main = "total ", cex.main = 1))

