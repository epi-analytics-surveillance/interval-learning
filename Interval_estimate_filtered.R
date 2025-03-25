####################################################################################
#	Script-file:   Interval_estimate_filtered.R
# Aim:           Using epidist packages to infer generation interval
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

# Set seed for random simulation 
set.seed(123)

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

# Estimate generation interval from filtered data
interval_estimate_filtered <- function(day=10){
  # Input: day, to only involve the data in which the infector infected times are before the input day
  
  # Converting the data to an epidist_linelist_data object
  data_linelist <- as_epidist_linelist_data(
    data$date_infector_infected, # lower bounds for primary times, infector infected day
    ptime_upr = data$date_infector_infected + 1, # upper bounds for primary times, infector infected day + 1
    stime_lwr = data$date_infectee_infected, # lower bounds for secondary times, infectee infected day
    stime_upr = data$date_infectee_infected + 1, # upper bounds for secondary times, infectee infected day + 1
    obs_time = max(data$date_infectee_infected) + 1 # max time in the dataset
  )
  
  # Prepare the data for use with the marginal model
  data_obs_prep <- as_epidist_naive_model(data_linelist) # using naive model instead of marginal model
  
  # Filter the data
  data_filtered <- filter(data_obs_prep, ptime_lwr <= day)
  
  # Fit the intercept only model
  intercept_only <- epidist(
    data = data_filtered,
    formula = mu ~ 1,
    family = exponential(), # choose exponential distribution
    algorithm = "sampling",
    chains = 2,
    cores = 2,
    refresh = ifelse(interactive(), 250, 0),
    seed = 1,
    backend = "cmdstanr"
  )
  
  # Plot posterior expectation of the delay distribution
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
  
  expectation_draws <- data_filtered |>
    data_grid(NA) |>
    add_marginal_dummy_vars() |>
    add_epred_draws(intercept_only, dpar = TRUE)
  
  epred_base_figure <- expectation_draws |>
    ggplot(aes(x = .epred)) +
    stat_halfeye() +
    labs(x = "", y = "", title = "Posterior Expectation of the Delay Distribution") +
    theme_minimal()
  
  # Plot continuous probability density function
  add_marginal_pdf_vars <- function(data) {
    return(
      mutate(
        data,
        relative_obs_time = Inf,
        pwindow = 0,
        swindow = 0,
        delay_upr = NA
      )
    )
  }
  
  draws_pdf <- data_filtered |>
    add_marginal_pdf_vars() |>
    add_predicted_draws(intercept_only, ndraws = 10)
  
  pdf_base_figure <- ggplot(draws_pdf, aes(x = .prediction)) +
    geom_density() +
    labs(x = "", y = "", title = "probability density function") +
    scale_x_continuous(limits = c(0, 30)) +
    theme_minimal()
  
  # Combine the figures
  print(epred_base_figure / pdf_base_figure)
  
 return(summary(intercept_only))
}
