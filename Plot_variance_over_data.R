####################################################################################
#	Script-file:   Plot_variance_over_data.R
# Aim:           To plot forward distribution: mean, sd, and variance, over data & date
# Note:          Filter the data a particular date, when both primary and secondary occurred before the date
# Project:       Capstone Project 2025
#
#	Data used:	   epidemic.csv
#	Data created:  none
# Student:       Zimo
####################################################################################

# Clean the session
rm(list=ls()) 

# Change working dictionary
setwd("~/Documents/Biostatistics/project/R code")

# Load the packages
library(haven)
library(epidist)
library(brms)
library(dplyr)
library(ggplot2)
library(tidybayes) 
library(modelr) 
library(cmdstanr)
library(patchwork)

# read in data
data <- read.csv("epidemic.csv")
summary(data)

set.seed(900915)

# Estimate generation interval from filtered data
mean_estimate_truncation <- function(percentage=1, day=80) {
  # Input: Percentage: range from 0 to 1; day: range from 1 to 80 
  # Percentage: to randomly filter down the data to set percentage (e.g., percentage = 0.1 will leave 10% data to use)
  # Day: to only involve the observations whose infector infected date are before the input day
  # Default value: day = 80 (i.e., maximum date), percentage = 1, complete data
  # If both arguments are set, the function will filter the data down firstly by percentage, then by date
  
  # Filter data down by percentage
  data_filtered <- data[sample(nrow(data), size = nrow(data) * percentage), ]
  
  # Filter the data down by date
  # data_filtered <- filter(data_filtered, date_infector_infected <= day)
  data_filtered <- filter(data_filtered, 
                          date_infector_infected <= day,
                          date_infectee_infected <= day)
  
  # Convert the data to an epidist_linelist_data object
  data_linelist <- as_epidist_linelist_data(
    data_filtered$date_infector_infected, # lower bounds for primary times, infector infected day
    ptime_upr = data_filtered$date_infector_infected + 1, # upper bounds for primary times, infector infected day + 1
    stime_lwr = data_filtered$date_infectee_infected, # lower bounds for secondary times, infectee infected day
    stime_upr = data_filtered$date_infectee_infected + 1, # upper bounds for secondary times, infectee infected day + 1
    obs_time = max(data_filtered$date_infectee_infected) + 1 # max time in the dataset
  )
  
  # Prepare the data for use with the naive model
  data_obs_prep <- as_epidist_naive_model(data_linelist) # using naive model instead of marginal model
  
  # Fit the intercept only model
  intercept_only <- epidist(
    data = data_obs_prep,
    formula = mu ~ 1,
    family = exponential(), # choose exponential distribution
    algorithm = "sampling",
    chains = 2,
    cores = 2,
    refresh = ifelse(interactive(), 250, 0),
    seed = 1,
    backend = "cmdstanr"
  )
  
  # Generate posterior distribution
  add_marginal_dummy_vars <- function(data) {
    return(
      mutate(
        data,
        relative_obs_time = NA,
        pwindow = NA,
        delay_upr = NA,
        swindow = NA
      )
    )
  }
  
  expectation_draws <- data_obs_prep |>
    data_grid(NA) |>
    add_marginal_dummy_vars() |>
    add_epred_draws(intercept_only, dpar = TRUE)
  
  g_mean <- mean(expectation_draws$.epred)
  g_sd <- sd(expectation_draws$.epred)
  g_ci <- quantile(expectation_draws$.epred, probs = c(0.025, 0.975))
  
  results <- list(
    obs = nrow(data_filtered),
    mean = g_mean,
    sd = g_sd,
    ci = g_ci
  )

  return(results)
}

###########################################################################################
# generate generation interval distribution plot by percentage
interval_plot_p <- function(day=80, step=0.1) {
  # Pre-allocate a posterior data frame
  posterior_df <- data.frame(i = 1:100, Mean = NA, SD = NA, lower_ci = NA, upper_ci = NA)
  
  # Define a sequence: 1, 5, 10, ..., 80
  i_values <- seq(from = 0.1, to = 1, by = step)
  
  # Generate the data frame that contains posterior distribution parameters for chosen days
  for (i in i_values) {
    result <- mean_estimate_truncation(i, day)
    # Assign values by index
    posterior_df$Mean[i*100] <- result$mean
    posterior_df$sd[i*100] <- result$sd
    posterior_df$var[i*100] <- (result$sd)^2
    posterior_df$lower_ci[i*100] <- result$ci["2.5%"]
    posterior_df$upper_ci[i*100] <- result$ci["97.5%"]
  }
  
  # Filter out the missing value and save the posterior data in the environment
  posterior_df <<- posterior_df[complete.cases(posterior_df), ]
  
  # Plot the mean generation interval with its standard deviation
  mean_interval_plot <- ggplot(posterior_df, aes(x = i, y = Mean)) +
    # Mean points
    geom_point() +  
    
    # Standard deviation (error bars)
    geom_errorbar(aes(ymin = Mean - sd, ymax = Mean + sd)) +
    
    # Labels and theme
    labs(x = "Data", y = "Mean Generation Interval (days)", 
         title = "Mean Generation Interval") +
    theme_minimal()
  
  # Plot the standard deviation 
  sd_plot <- ggplot(posterior_df, aes(x = i, y = sd)) +
    # Sd points
    geom_point() + 
    
    # Smooth line
    geom_smooth(method = "loess") +
    
    # Labels and theme
    labs(x = "Percentage of Data (%)", y = "Generation Interval Standard Deviation", 
         title = "Generation Interval Standard Deviation") +
    theme_minimal()
  
  # Plot the variance 
  var_plot <- ggplot(posterior_df, aes(x = i, y = var)) +
    # Variance points
    geom_point() +  
    
    # Smooth line
    geom_smooth(method = "loess") +
    
    # Labels and theme
    labs(x = "Percentage of Data (%)", y = "Generation Interval Variance", 
         title = "Generation Interval Variance") +
    theme_minimal()
  
  # Combine the result plot  
  result_plot <- list(mean = mean_interval_plot,
                      sd = sd_plot,
                      var = var_plot)
  return(result_plot)
}

interval_plot_p(80,0.01)

#################################################################################
# Generate forward distribution over time
interval_plot <- function(percentage=1, day=80) {
  # Pre-allocate a posterior data frame
  posterior_df <- data.frame(i = 1:80, Obs = NA, Mean = NA, SD = NA, lower_ci = NA, upper_ci = NA)
  
  # Define a sequence: 1, 5, 10, ..., 80
  i_values <- seq(from = 1, to = day, by = 1)
  
  # Generate the data frame that contains posterior distribution parameters for chosen days
  for (i in i_values) {
    result <- mean_estimate_truncation(percentage, i)
    # Assign values by index
    posterior_df$Obs[i] <- result$obs
    posterior_df$Mean[i] <- result$mean
    posterior_df$SD[i] <- result$sd
    posterior_df$lower_ci[i] <- result$ci["2.5%"]
    posterior_df$upper_ci[i] <- result$ci["97.5%"]
  }
  
  # Filter out the missing value
  plot_data <<- posterior_df[complete.cases(posterior_df), ]
  
  # Plot the mean generation interval with its standard deviation
  mean_interval_plot <- ggplot(plot_data, aes(x = i, y = Mean)) +
    # Mean points
    geom_point() +  
    
    # Standard deviation (error bars)
    geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD)) +
    
    # Labels and theme
    labs(x = "Time (days)", y = "Mean Generation Interval (days)", 
         title = "Mean Generation Interval Over Time (Days)") +
    theme_minimal()
  
  return(mean_interval_plot)
}

interval_plot()

