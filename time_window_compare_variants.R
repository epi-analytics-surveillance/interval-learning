####################################################################################
#	Script-file:   time_window_compare_variants.R
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

# Read in data
data <- read.csv("combined_epidemic.csv")
summary(data)

################################################################################
# Initial exploration
data$generation_interval <- data$date_infectee_infected - data$date_infector_infected
data$pathogen <- as.factor(data$pathogen)

# Function to group, summarise, and plot the count, by the infector infected time
count_plot <- function(data=data, time_interval = 10){
  # Input: time-interval: Group data into time_intervals bins 
  # Initial value is timer_interval = 10, meaning days 1–10 = Group 1, 11–20 = Group 2, ..., 71–80 = Group 8
  
  # Group data by input interval
  data$time_group <- cut(data$date_infector_infected, 
                         breaks = seq(min(data$date_infector_infected), 
                                      max(data$date_infector_infected)+1, 
                                      by = time_interval),
                         right = FALSE, 
                         include.lowest = TRUE)
  
  # Generate summary statistics of generation interval for each group
  summary_stats <- data %>%
    group_by(time_group,pathogen) %>%
    summarise(
      count = n(),
      mean = mean(generation_interval),
      median = median(generation_interval),
      sd = sd(generation_interval),
      variance = var(generation_interval),
      q25 = quantile(generation_interval, 0.25),
      q75 = quantile(generation_interval, 0.75),
      min = min(generation_interval),
      max = max(generation_interval)
    )
  print(summary_stats, n = Inf)  # Print all rows
  
  # Count plot
  plot <- ggplot(summary_stats, aes(x = time_group, y = count, color = pathogen, group = pathogen)) +
    geom_smooth(method = "loess", se = F, linewidth = 1.2, span = 0.7) +
    geom_point(size = 3) + # Show original data points
    labs(x = "Time (days)", y = "Count", fill = "Pathogen",
      title = paste("Number of Infectors per", time_interval, "Days by Pathogen")) +
    theme_minimal()
  
  print(plot)
}

count_plot(data,10)

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
time_window_plot <- function(data=data, percentage=1, step=1, time_window=10, secondary_interval=40,
                             start_day, end_day) {
  #Input: step: the distance between two time window;
  # E.g., default = 1, meaning the start day of the next time window is 1 day after the start day of the previous one
  
  # Pre-allocate a posterior data frame 
  posterior_df <- data.frame(i = seq(0, end_day-start_day, 1), Obs = NA,
                             primary_day = NA, primary_end = NA, secondary_end = NA,
                             Mean = NA, SD = NA, lower_ci = NA, upper_ci = NA)
  
  
  # Define a sequence for the primary_day
  i_values <- seq(from = start_day, to = end_day, by = step)
  
  # Generate the data frame that contains posterior distribution parameters for chosen days
  for (i in i_values) {
    # Estimate the parameters from the previous function
    result <- mean_estimate_time_window(data=data, percentage=percentage, 
                                        primary_day=i, time_window=time_window,
                                        secondary_interval=secondary_interval)
    
    # Check if filtered data is empty
    if (is.atomic(result)) {
      posterior_df[i-start_day+1, ] <- NA
      next
    }
    
    # Extract results and write in to the data frame
    posterior_df$Obs[i-start_day+1] <- result$obs
    posterior_df$primary_day[i-start_day+1] <- result$p_day
    posterior_df$primary_end[i-start_day+1] <- result$p_end
    posterior_df$secondary_end[i-start_day+1] <- result$s_end
    posterior_df$Mean[i-start_day+1] <- result$mean
    posterior_df$SD[i-start_day+1] <- result$sd
    posterior_df$lower_ci[i-start_day+1] <- result$ci["2.5%"]
    posterior_df$upper_ci[i-start_day+1] <- result$ci["97.5%"]
  }
  
  # Select the complete data
  posterior_df <<- posterior_df[complete.cases(posterior_df),]
  
  dataset_name <- deparse(substitute(data)) # Extract data name
  
  # Plot the mean generation interval with its standard deviation
  plot <- ggplot(posterior_df, aes(x = primary_day, y = Mean)) +
    # Mean points
    geom_point() +  
    
    # Standard deviation (error bars)
    geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD)) +
    
    # Labels and theme
    labs(x = "Start Date of Time Window", y = "Mean Generation Interval (days)", 
         title = paste("Mean Generation Interval using", percentage*100, "% data for", dataset_name),
         caption = paste("Note: Each time window is", time_window, "days")) +
    
    theme_minimal()
  
  return(plot)
}

v0 <- time_window_plot(data, 1, 1, 10, 40, 0, 90)

# Plot the mean generation interval separately for pathogen 1 and pathogen 2 
pathogen_1 <- data[data$pathogen==1,]
v1 <- time_window_plot(pathogen_1, 1, 1, 10, 40, 0, 60)
pathogen_2 <- data[data$pathogen==2,]
v2 <- time_window_plot(pathogen_2, 1, 1, 10, 40, 20, 90)

# Combine the plots
combined_plot <- v0 / v1 / v2 +
  plot_annotation(title = "Generation Interval Distribution for complete data, pathogen 1 and pathogen 2", 
                  theme = theme(plot.title = element_text(hjust = 0.5)))

################################################################################
# Combine the result data set for pathogen 1 and 2
combine_data_pathogen <- function(data=data, p1_start, p1_end, p2_start, p2_end){
  
  # Compute result for pathogen 1
  pathogen_1 <- data[data$pathogen==1,]
  p1 <- time_window_plot(pathogen_1, 1, 1, 10, 40, p1_start, p1_end)
  data1 <- posterior_df
  data1$pathogen <- "1"
  
  # Compute result for pathogen 2
  pathogen_2 <- data[data$pathogen==2,]
  p2 <- time_window_plot(pathogen_2, 1, 1, 10, 40, p2_start, p2_end)
  data2 <- posterior_df
  data2$pathogen <- "2" 
  
  # Combine the results together
  combined_result <- bind_rows(data1, data2)
  
  return(combined_result)
}

# Time window plot for pathogen1 from day 0 to 10, pathogen2 from day 20 to 30
data_plot1 <- combine_data_pathogen(data, 0, 10, 20, 30)

ggplot(data_plot1, aes(x = i, y = Mean)) +
  geom_point() +  # Mean points
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD)) + # Standard deviation (error bars)
  labs(x = "Relative Time", y = "Mean Generation Interval (days)", 
       title = "Initial 10 Days for Pathogen 1 (day 0-10) and 2 (day 20-30)",
       caption = "Note: Each time window is 10 days") +
  scale_x_continuous(breaks = seq(0,10,1)) +
  facet_wrap(~pathogen) +
  theme_minimal()


# Time window plot for pathogen1 from day 20 to 30, pathogen2 from day 20 to 30
data_plot2 <- combine_data_pathogen(data, 20, 30, 20, 30)

ggplot(data_plot2, aes(x = primary_day, y = Mean)) +
  geom_point() +  # Mean points
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD)) + # Standard deviation (error bars)
  scale_x_continuous(breaks = seq(20,30,1)) +
  labs(x = "Start Date of Time Window", y = "Mean Generation Interval (days)", 
       title = "Pathogen 1 and 2 at Absolute Time Day 20-30",
       caption = "Note: Each time window is 10 days") +
  facet_wrap(~pathogen) +
  theme_minimal()


