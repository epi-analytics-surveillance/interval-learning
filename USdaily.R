################################################################################
#	Script-file:   USdaily.R
# Aim:           Contour Plot
# Project:       Capstone Project 2025
#
#	Data used:	   US.csv
#	Data created:  data
# Student:       Zimo
################################################################################

# Clean the session
rm(list=ls()) 

# Change working dictionary
setwd("~/Documents/Biostatistics/project/R code")

library(tidyverse)
library(ggplot2)
library(patchwork)
library(haven)

realdata <- read_csv("US.csv")

# Select the data from 2020.01.23 to 2021.08.19 as used in the paper
realdata$date <- as.Date(realdata$date)
realdata <- realdata[realdata$date >= as.Date("2020-01-23") & realdata$date <= as.Date("2021-08-25"), ]
realdata <- realdata[realdata$overall_outcome == "Positive",]

# Plot figure.2 in the paper
incidence_plot <- list()

for (i in c("Alabama", "Arizona", "Iowa", "Michigan", "New Hampshire", "New Mexico")) {
  state <- realdata[realdata$state_name == i,]
  
  p <- ggplot(state, aes(date, new_results_reported)) +
    geom_line() +
    labs(title = i,
         x = "Date",
         y = "Cases per Million People",
         caption = " CDC statewide incidence data from 1/23/20 to 08/25/21") + 
    theme_bw()

  incidence_plot[[i]] <- p
}

# Combine the plots together
incidence_plot <- discard(incidence_plot, is.null)
combined_incidence_plot <- wrap_plots(incidence_plot, ncol = 3) # 3 columns
print(combined_incidence_plot)

# case count
state <- realdata[realdata$state_name == "Alabama",]
case_count <- state$new_results_reported[state$overall_outcome == "Positive"]

################################################################################
# Log-likelihood function
log_likelihood <- function(I, R_t, W_s) {
  # I: Daily observed case counts (length N)
  # R_t: Reproduction numbers (length N)
  # W_s: Generation interval probabilities (length <= N)
  
  n_days <- length(I) # Number of days == length(R_t)
  E_I <- numeric(n_days)
  E_I[1] <- I[1] # Start with I[1] cases
  ll <- I[1] * log(E_I[1]) - E_I[1] - lgamma(I[1] + 1) # Initial value for log-likelihood
  
  # Calculate expected cases lambda_t
  for (t in 2:n_days) {
    
    E_I[t] <- 0
    for (s in 1:(t-1)) {
      
      # End if exceeds the length of W_s 
      if (s > length(W_s)) {
        break
      }
      
      E_I[t] <- E_I[t] + R_t[t - s] * W_s[s] * I[t - s]
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
  meanlog <- log((GI_mean^2) / sqrt(GI_mean^2 + GI_sd^2))
  
  # Gamma parameter for reproduction number
  shape <- (Rt_mean^2) / (Rt_sd^2)
  scale <- (Rt_sd^2) / Rt_mean
  
  # Generate GI distribution
  W_s_raw <- dlnorm(1:50, meanlog, sdlog)
  W_s <- W_s_raw/sum(W_s_raw) # standardise to make the sum equal to 1
  W_s_raw <- W_s_raw[W_s_raw  > 1e-5]
  
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
# Generate the whole data set function
ll_plot_data <- function(days, I,
                         Rt_lower, Rt_upper, Rt_step, GI_sd,
                         GI_lower, GI_upper, GI_step, Rt_sd){
  # N: the number of observations of mean, sd, and log-likelihood
  # case_count_days: the length of the simulated case count data
  # I: Simulated series of case count
  # Rt & GI: from lower to upper, by step
  
  # generate a data frame
  GI_mean <- seq(GI_lower,GI_upper, GI_step) 
  Rt_mean <- seq(Rt_lower,Rt_upper, Rt_step)
  df <- expand.grid(GI_mean = GI_mean, Rt_mean = Rt_mean)
  N <- nrow(df) # number of rows == length(GI_mean)*length(Rt_mean)
  df$GI_sd = rep(GI_sd, N)
  df$Rt_sd = rep(Rt_sd, N)
  df$Rt_length = rep(days, N)
  
  # generate the GI and Rt based on the simulated mean and sd
  results <- apply(df, 1, function(row) {
    do.call(sim_GI_Rt, as.list(row)) # apply the function sim_GI_Rt to each row
  })
  
  # add the generated GI and Rt to the data set 
  results_df <- map_dfr(results, as_tibble)
  final_data <- bind_cols(df, results_df)
  
  # Calculate the log-likelihood for each row
  for (i in 1:N) {
    final_data$ll[i] <- log_likelihood(I, final_data$Rt[[i]], final_data$Ws[[i]])
  }
  
  return(final_data)
}

################################################################################
# Set parameters
parameters <- data.frame(
  GI_mean = 7,
  GI_sd = 1,
  GI_lower = 4.5,
  GI_upper = 8,
  Rt_mean = 1,
  Rt_sd = 0.01,
  Rt_lower = 0.8,
  Rt_upper = 1.3,
  days = length(case_count)
)

# Calculate log-likelihood
data <- ll_plot_data(days = parameters$days, I = case_count,
                     Rt_lower = parameters$Rt_lower, Rt_upper = parameters$Rt_upper, 
                     Rt_step = 0.1, Rt_sd = parameters$Rt_sd,
                     GI_lower = parameters$GI_lower, GI_upper = parameters$GI_upper, 
                     GI_step = 0.1, GI_sd = parameters$GI_sd)

data <- data[data$ll  >= -200000,] # zoom in 

# Generate contour plot
ggplot(data, aes(x = GI_mean, y = Rt_mean, z = ll)) +
  geom_contour_filled(bins=10) +  # Filled contours
  geom_contour() +  # Add contour lines
  scale_fill_viridis_d(name = "Log-Likelihood") +   # Color scale
  labs(title = "Log-Likelihood Contour Plot",
       x = "Generation Interval Mean (days)",
       y = "Reproduction Number",
       caption = paste0("State: ", "Alaska")) + # case count as caption
  theme(plot.caption = element_text(hjust = 0)) + # location of caption
  scale_x_continuous(limits = c(4.5, 8)) + # x scale range
  #scale_y_continuous(limits = c(0.8, 1.2)) + # y scale range
  geom_point(data=parameters, aes(x = GI_mean, y = Rt_mean), # mark true value
             size = 5, shape = 10, color = "black", inherit.aes = FALSE) +
  geom_text(data=parameters, aes(x = GI_mean, y = Rt_mean), # label true value
            label=paste0("(",parameters$GI_mean,",",parameters$Rt_mean,")"), 
            nudge_y = 0.03, color = "black", inherit.aes = FALSE) +
  theme_minimal()

