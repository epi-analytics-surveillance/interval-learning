####################################################################################
#	Script-file:   Plot_mean_generation_interval.R
# Aim:           To estimate generation interval parameters and plot its mean with sd
# Project:       Capstone Project 2025
#
#	Data used:	   epidemic.csv
#	Data created:  none
# Student:       Zimo
####################################################################################

# Clean the session
#rm(list=ls()) 

# Change working dictionary
#setwd("~/Documents/Biostatistics/project/R code")

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
interval_mean_estimate <- function(percentage=1, day=80) {
  # Input: Percentage: range from 0 to 1; day: range from 1 to 80 
  # Percentage: to randomly filter down the data to set percentage (e.g., percentage = 0.1 will leave 10% data to use)
  # Day: to only involve the observations whose infector infected date are before the input day
  # Default value: day = 80 (i.e., maximum date), percentage = 1, complete data
  # If both arguments are set, the function will filter the data down firstly by percentage, then by date
  
  # Filter data down by percentage
  data_filtered <- data[sample(nrow(data), size = nrow(data) * percentage), ]
  
  # Filter the data down by date
  data_filtered <- filter(data_filtered, date_infector_infected <= day)
  
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

# Generate mean generation interval plot
interval_plot <- function(percentage=1, day=80, step=1) {
  # Pre-allocate a posterior data frame
  posterior_df <- data.frame(i = 1:80, Obs = NA, Mean = NA, SD = NA, lower_ci = NA, upper_ci = NA)
 
   # Define a sequence: 1, 5, 10, ..., 80
  i_values <- seq(from = 1, to = day, by = step)
  
  # Generate the data frame that contains posterior distribution parameters for chosen days
  for (i in i_values) {
    result <- interval_mean_estimate(percentage, i)
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


# Initialize lists to store plots and data
plot_list <- list()
data_list <- list()

# Generate plots using complete, 75%, 50%, and 25% data for the first 20 days
for (i in c(1, 0.75, 0.5, 0.25)) {
  # Generate plot and data
  p <- interval_plot(i, 20, 1)
  plot_data <- plot_data  # Assuming interval_plot stores data in plot_data
  
  # Store results with descriptive names
  plot_name <- paste0("p_", i * 100)
  data_name <- paste0("data_", i * 100)
  
  plot_list[[plot_name]] <- p
  data_list[[data_name]] <- plot_data
}

# Combine the plots together
plot_list$p_100 / plot_list$p_75 / plot_list$p_50 / plot_list$p_25
