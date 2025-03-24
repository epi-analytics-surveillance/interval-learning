####################################################################################
#	Script-file:   Explore_simulated_epidemic.R
# Aim:           Extract summarised statistics from simulated case tracing data
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
library(ggplot2)
library(dplyr)
library(patchwork)

# Read in the simulated data
data <- read.csv("epidemic.csv")

# Explore the data
summary(data)

# Calculate generation interval
data$generation_interval <- data$date_infectee_infected - data$date_infector_infected

# Plot the generation interval over time to explore trend
plot_GI <- ggplot(data, aes(x = date_infector_infected, y = generation_interval)) +
  geom_point(alpha = 0.9) +  # Scatter plot
  labs(title = "Scatter Plot of Generation Interval over Time", x = "Time (days)", y = "Generation Interval (days)") +
  theme_minimal()
plot_GI # Print the plot

# Categorize the data into groups based on infectors' infection onset time
data$time_group <- cut(data$date_infector_infected, 
    breaks = seq(0, 80, by = 10),
    right = FALSE, 
    include.lowest = TRUE)


# Function to group, summarise, and plot data
summary_stats <- function(time_interval = 10){
  # Input: time-interval: Group data into time_intervals bins 
  # Initial value is timer_interval = 10, meaning days 1–10 = Group 1, 11–20 = Group 2, ..., 71–80 = Group 8
  
  # Group data by input interval
  data$time_group <- cut(data$date_infector_infected, 
                         breaks = seq(0, 80, by = time_interval),
                         right = FALSE, 
                         include.lowest = TRUE)
  
  # Generate summary statistics of generation interval for each group
  summary_stats <- data %>%
    group_by(time_group) %>%
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
  
  # Plot mean of generation intervals across group 
  plot_meansd <- ggplot(summary_stats, aes(x = time_group, y = mean, group = 1)) +
    geom_point() +  # Scatter plot
    geom_line() +
    geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.2) +
    labs(title = "Trend of Mean ± SD Across Groups", x = "Time (days)", y = "Mean Generation Interval (days)") +
    theme_bw()
  
  # Plot median with 25 and 75 percentile
  plot_median <- ggplot(summary_stats, aes(x = time_group, y = median, group = 1)) +
    geom_point() +  # Scatter plot
    geom_line() +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.2) +
    labs(title = "Trend of Median with 25-75 percentile Across Groups", x = "Time (days)", y = "Median Generation Interval (days)") +
    theme_bw()
  
  # Combine plots using patchwork
  combined_plot <- plot_meansd / plot_median  # Stack plots vertically
  print(combined_plot)
}
