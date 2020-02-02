refit_p_confirm_delay <- TRUE # if TRUE, fit geometric distribution to confirmation delay data;
# if FALSE, read from file
bayesian_p_confirm_delay <- FALSE # if TRUE, use posterior for confirmation delay parameter, if FALSE, use point estimate

#setwd("~/Documents/case_to_infection/")
# setwd("~/GitHub/case_to_infection/")

library(ggplot2)
library(tidyverse)
library(lubridate)
library(cowplot)
library(patchwork)
library(ggpmisc)
library(ggpubr)
library(maptools)
library(maps)
library(data.table)
library(googlesheets4)
if(refit_p_confirm_delay && bayesian_p_confirm_delay) {
  library(rstan)
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)
}
source("code/analysis_functions.R")

symptom_times <- read.csv("augmented_data/inflated_symptom_times.csv")
symptom_times <- symptom_times[,1] %>% unlist %>% as.Date(origin = "1970-01-01")

# date_today <- convert_date(Sys.Date())
date_today <- as.Date("2020-02-01") # because Ada last pulled the data on 2020-02-01 -- to fix

weibull_stan_draws <- read.csv("data/backer_weibull_draws.csv")

## source("code/plot_china_map.R")

## These column names will be kept as keys for the rest of the analysis
key_colnames <- c("ID","age","country", "sex","city","province",
                  "latitude","longitude")
## These are the column names with used variables
var_colnames <- c("date_confirmation","date_onset_symptoms","date_admission_hospital")
use_colnames <- c(key_colnames, var_colnames)

## Number of bootstrap samples to take. Set this to something small for a quick run
repeats <- 1000

## load the data - try to only do this once otherwise auth token gets stale
## First step is to clean and take a look at the data
## This combines the data for Hubei and other locations in China
## Note this is ONLY China

# kudos_dat <- readRDS("data/kudos_dat.rds")

# check if we have already loaded the objects, load if not
if(!exists("kudos_dat")) {
  source("code/pull_kudos_linelist.R")
}

# saveRDS(kudos_dat, "data/kudos_dat.rds")
####################################
## KUDOS LINE LIST CONFIRMATION DELAY
####################################
## Fit a geometric distribution to the confirmation delay distribution
if(refit_p_confirm_delay) {
  use_delays_kudos <- kudos_dat %>% select(delay) %>% drop_na() %>% pull(delay)
  if(bayesian_p_confirm_delay) { # bayesian fit (posterior)
    confirmation_delay_model <- stan_model("code/geometric.stan")
    fit_kudos <- fit_geometric_stan(use_delays_kudos - 1, confirmation_delay_model) %>%
      extract(pars = "p")
    names(fit_kudos) <- "par"
  } else {
    # frequentist fit (point estimate)
    fit_kudos <- optim(c(0.1), fit_geometric, dat=use_delays_kudos-1,method="Brent",lower=0,upper=1)
  }
} else {
  ## or read from file
  if(bayesian_p_confirm_delay) {
    fit_kudos <- read.csv("data/p_confirm_delay_draws.csv") %>%
      as.list
    names(fit_kudos) <- "par"
  } else {
    fit_kudos <- list(par = as.numeric(read.csv("data/p_confirm_delay.csv")))
  }
}

plot_times <- seq(0,max(kudos_dat$delay,na.rm=TRUE))
predict_delay <- function(p_confirm_delay) {
  dgeom(plot_times, prob = p_confirm_delay)
} 

p_confirm_delay_kudos <- kudos_dat %>% select(delay) %>% drop_na() %>%
  ggplot() + 
  geom_histogram(aes(x=delay,y=..density..),binwidth=1,col="black") +
  scale_x_continuous(breaks=seq(0,max(kudos_dat$delay,na.rm=TRUE),by=5),labels=seq(0,max(kudos_dat$delay,na.rm=TRUE),by=5)) +
  scale_y_continuous(expand=c(0,0),limits=c(0,0.15)) +
  geom_vline(xintercept=1,linetype="dashed") +
  ylab("Probability density") + xlab("Days since symptom onset") +
  ggtitle("Distribution of delays between symptom onset\n and confirmation, Kudos line list data") +
  theme_pubr()
if(bayesian_p_confirm_delay) {
  fit_kudos_line <- vapply(fit_kudos$par, predict_delay, numeric(length(plot_times))) %>%
    apply(1, quantile, probs = c(0.025, 0.975))
  fit_line_kudos_dat <- data.frame(x=plot_times + 1,ymin=fit_kudos_line[1,],
                                   ymax = fit_kudos_line[2,])
  
  p_confirm_delay_kudos <- p_confirm_delay_kudos +
    geom_ribbon(data = fit_line_kudos_dat, aes(x = x, ymin = ymin, ymax = ymax), fill = "red")
} else {
  fit_kudos_line <- dgeom(seq(0,max(kudos_dat$delay,na.rm=TRUE),by=1),prob=fit_kudos$par)
  fit_line_kudos_dat <- data.frame(x=seq(1,max(kudos_dat$delay,na.rm=TRUE)+1,by=1),y=fit_kudos_line)
  
  p_confirm_delay_kudos <- p_confirm_delay_kudos +
    geom_line(data=fit_line_kudos_dat, aes(x=x,y=y), col="red",size=1)
}

p_confirm_delay_kudos
sim_data_confirmation <- matrix(NA, nrow=repeats, ncol=length(symptom_times))

for(i in seq_len(repeats)){
  # sample from posterior if bayesian
  if(bayesian_p_confirm_delay) {
    p_confirm_delay <- sample(fit_kudos$par,1)
  } else {
    # use point estimate if frequentist
    p_confirm_delay <- fit_kudos$par
  }
  
  ## Get symptom onset and confirmation times
  
  sim_data_confirmation[i,] <- simulate_confirmation_times(symptom_times, 
                                                              p_confirm_delay=p_confirm_delay)
}

sim_data_confirmation_melted <- reshape2::melt(sim_data_confirmation)

sim_data_confirmation_melted$var <- "date_confirmation"

colnames(sim_data_confirmation_melted) <- c("repeat_no","individual","date","var")

##  conert to dates
sim_data_confirmation_melted$date <- as.Date(floor(sim_data_confirmation_melted$date), origin="1970-01-01")

## Sum by repeat, variable and date ie. events per day
sim_data_sum <- sim_data_confirmation_melted %>% group_by(repeat_no, var, date) %>% tally()
sim_data_sum <- sim_data_sum %>% group_by(repeat_no, var) %>%
  mutate(date_diff = as.numeric(date_today - date))

## OVERALL PLOT
## Distribution of times for each date
sim_data_quantiles <- sim_data_sum %>% group_by(date, var) %>% 
  do(data.frame(t(c(quantile(.$n, probs = c(0.025,0.5,0.975),na.rm=TRUE),mean(.$n))))) %>% 
  arrange(date)

colnames(sim_data_quantiles) <- c("date","Variable","lower","median","upper","mean")
onset_data <- data.frame(date_onset_symptoms = symptom_times) %>%
  group_by(date_onset_symptoms) %>% tally()
# onset_data <- combined_dat_final %>% filter(!is.na(date_onset_symptoms)) %>% group_by(date_onset_symptoms) %>% tally()
# onset_data$Variable <- "onsets"

sim_data_quantiles_truncated <- sim_data_quantiles %>% filter(date <= convert_date(Sys.Date()))
sim_data_plot <- plot_forward_simulation(sim_data_quantiles_truncated, onset_data,ymax=5000,ybreaks=500,
                                           max_date = date_today, min_date="01.01.2020")
sim_data_plot