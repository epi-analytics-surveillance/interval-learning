####################################################################################
#	Script-file:   Infer_generation_interval.R
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
library(tidybayes) 
library(modelr) 
library(cmdstanr)

# read in data
data <- read.csv("epidemic.csv")

# Converting the data to an epidist_linelist_data object
data_linelist <- as_epidist_linelist_data(
  data$date_infector_infected, # lower bounds for primary times, infector infected day
  ptime_upr = data$date_infector_infected + 1, # upper bounds for primary times, infector infected day + 1
  stime_lwr = data$date_infectee_infected, # lower bounds for secondary times, infectee infected day
  stime_upr = data$date_infectee_infected + 1, # upper bounds for secondary times, infectee infected day + 1
  obs_time = max(data$date_infectee_infected) + 1 # max time in the dataset
)

# Prepare the data for use with the marginal model
data_obs_prep <- as_epidist_marginal_model(data_linelist)

# Fit the intercept only model
fit <- epidist(
  data = data_obs_prep,
  formula = mu ~ 1,
  family = lognormal(),
  algorithm = "sampling",
  chains = 2,
  cores = 2,
  refresh = ifelse(interactive(), 250, 0),
  seed = 1,
  backend = "cmdstanr"
)
summary(fit)
