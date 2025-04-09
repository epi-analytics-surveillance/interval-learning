####################################################################################
#	Script-file:   time_window.R
# Aim:           To plot generation interval distribution with in time windows
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
data <- read.csv("epidemic2.csv")
summary(data)

################################################################################
mean_estimate_time_window <- function(data=data, percentage=1, primary_day=0, time_window=10, secondary_interval = 40) {
  # Input: Percentage: range from 0 to 1; primary_day: range from 0 to 75, 
  # secondary_interval range from 1 to 80, time_window range from 1 to 80,
  
  # percentage: to randomly filter down the data to set percentage (e.g., percentage = 0.1 will leave 10% data to use)
  # primary_day: to filter the data by time, with primary cases occurred after the primary_day
  # time_window: to filter the data by time, with primary cases occurred before the primary_day + time_window
  # secondary_interval: to filter the data by time, with secondary cases occurred before the primary_day + secondary_interval
  
  # Filter the data down firstly by percentage, then by date
  # Set the date when the latest primary cases occurred
  primary_max <- max(data$date_infector_infected)
  secondary_max <- max(data$date_infectee_infected)
  
  # Set the date that the time_interval ends to be smaller than the latest date when primary cases occurred
  primary_cut <- pmin(primary_day+time_window, primary_max)
  secondary_cut <- pmin(primary_day+secondary_interval, secondary_max)
  
  # Filter data down by percentage
  data_filtered <- data[sample(nrow(data), size = nrow(data) * percentage), ]
  
  # Filter the data down by date
  data_filtered <- filter(data_filtered, 
                          date_infector_infected >= primary_day,
                          date_infector_infected < primary_cut,
                          date_infectee_infected < secondary_cut)
  
  # Check if filtered data is empty
  if (nrow(data_filtered) == 0) {
    return(NA)
  }
  
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
    p_day = primary_day,
    p_end = primary_cut,
    s_end = secondary_cut,
    obs = nrow(data_filtered),
    mean = g_mean,
    sd = g_sd,
    ci = g_ci
  )

  return(results)
}

################################################################################
# Generate mean generation interval plot over time window
time_window_plot <- function(data=data, percentage=1, step=1, time_window=10, secondary_interval=40) {
  #Input: step: the distance between two time window;
  # E.g., default = 1, meaning the start day of the next time window is 1 day after the start day of the previous one
  
  # Set the date when the latest primary cases occurred
  primary_max <- max(data$date_infector_infected)
  
  # Pre-allocate a posterior data frame 
  posterior_df <- data.frame(i = 1:primary_max, Obs = NA,
                             primary_day = NA, primary_end = NA, secondary_end = NA,
                             Mean = NA, SD = NA, lower_ci = NA, upper_ci = NA)

  
  # Define a sequence for the primary_day
  i_values <- seq(from = 0, to = primary_max - time_window, by = step)
  
  # Generate the data frame that contains posterior distribution parameters for chosen days
  for (i in i_values) {
    # Estimate the parameters from the previous function
    result <- mean_estimate_time_window(data=data, percentage=percentage, 
                                        primary_day=i, time_window=time_window,
                                        secondary_interval=secondary_interval)
    
    # Check if filtered data is empty
    if (is.atomic(result)) {
      posterior_df[i + 1, ] <- NA
      next
    }
    
    # Extract results and write in to the data frame
    posterior_df$Obs[i+1] <- result$obs
    posterior_df$primary_day[i+1] <- result$p_day
    posterior_df$primary_end[i+1] <- result$p_end
    posterior_df$secondary_end[i+1] <- result$s_end
    posterior_df$Mean[i+1] <- result$mean
    posterior_df$SD[i+1] <- result$sd
    posterior_df$lower_ci[i+1] <- result$ci["2.5%"]
    posterior_df$upper_ci[i+1] <- result$ci["97.5%"]
  }
  
  # Select the complete data
  posterior_df <<- posterior_df[complete.cases(posterior_df),-which(names(posterior_df) == "i")]
    
  # Plot the mean generation interval with its standard deviation
  plot <- ggplot(posterior_df, aes(x = primary_day, y = Mean)) +
    # Mean points
    geom_point() +  
    
    # Standard deviation (error bars)
    geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD)) +
    
    # Labels and theme
    labs(x = "Start Date of Time Window", y = "Mean Generation Interval (days)", 
         title = paste("Mean Generation Interval using", percentage*100, "% data"),
         caption = paste("Note: Each time window is", time_window, "days")) +
    
    theme_minimal()
  
  return(plot)
}

# Initialize lists to store plots and data
plot_list <- list()
data_list <- list()

# Generate plots using complete, 1%, 10%, 25%, 50%, 75%, and complete data
for (i in c(0.01, 0.1, 0.25, 0.5, 0.75, 1)) {
  # Generate plot and data
  p <- time_window_plot(data, i, 1, 10, 40)
  plot_df <- posterior_df  # Assuming interval_plot stores data in plot_data
  
  # Store results with descriptive names
  plot_name <- paste0("p_", i * 100)
  data_name <- paste0("data_", i * 100)
  
  plot_list[[plot_name]] <- p          # Access the result: e.g., plot_list$p_50
  data_list[[data_name]] <- plot_df    # Access the result: e.g., data_list$data_50
}

# Combine the plots together
combined_plot <- (plot_list$p_1 + plot_list$p_10 + plot_list$p_25) / 
  (plot_list$p_50 + plot_list$p_75 + plot_list$p_100) +
  plot_annotation(title = "Generation Interval Distribution at Different Percentage of Data", 
                  theme = theme(plot.title = element_text(hjust = 0.5)))

