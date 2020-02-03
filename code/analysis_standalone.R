setwd("~/Documents/case_to_infection/")
#setwd("~/GitHub/case_to_infection/")
refit_p_confirm_delay <- TRUE # if TRUE, fit geometric distribution to confirmation delay data;
# if FALSE, read from file
bayesian_p_confirm_delay <- TRUE # if TRUE, use posterior for confirmation delay parameter, if FALSE, use point estimate

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
library(extraDistr)

## Some setup for STAN
if(refit_p_confirm_delay && bayesian_p_confirm_delay) {
  library(rstan)
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)
}
source("code/analysis_functions.R")


date_today <- convert_date(Sys.Date())

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
# combined_dat <- readRDS("data/combined_dat.rds")
# use_data_diff <- readRDS("data/use_data_diff.rds")

# check if we have already loaded the objects, load if not
if(!exists("combined_dat")) {
  source("code/pull_and_clean_linelist.R")
}
if(!exists("kudos_dat")) {
  source("code/pull_kudos_linelist.R")
}
if(!exists("use_data_diff")) {
  source("code/pull_arcgis_data.R")
}

# saveRDS(kudos_dat, "data/kudos_dat.rds")
# saveRDS(combined_dat, "data/combined_dat.rds")
# saveRDS(use_data_diff, "data/use_data_diff.rds")

## Key data objects produce:
## kudos_dat: the kudos line list data
## combined_dat: the Mortiz Kraemer line list data
## use_data_diff: the arcgis new confirmed case data

## The main thing returned from this is "combined_dat"
## Let's make combined_dat into confirmed case totals
## Date to expand over:
first_date <- min(combined_dat$date_confirmation,na.rm=TRUE)
last_date <- max(combined_dat$date_confirmation,na.rm=TRUE)
dates <- first_date:last_date
dates <- convert_date(dates)
all_combos <- expand.grid(province=as.character(unique(combined_dat$province)),
                          date_confirmation=dates) %>% arrange(province, date_confirmation)
all_combos$province <- as.character(all_combos$province)
confirmed_cases_linelist <- combined_dat %>% filter(!is.na(date_confirmation)) %>% 
  group_by(province, date_confirmation) %>% tally()  %>% ungroup() 

confirmed_cases_linelist <- confirmed_cases_linelist %>% right_join(all_combos) %>% 
  arrange(province, date_confirmation) %>% mutate(n=ifelse(is.na(n), 0, n))

confirmed_cases_linelist$province <- factor(confirmed_cases_linelist$province, levels=levels(combined_dat$province))
confirmed_cases_linelist <- confirmed_cases_linelist %>% 
  mutate(pre_reports=date_confirmation <= convert_date("21.01.2020"))

## Have a quick look to see how this compares to the ARCGIS data
p_confirm_pre_arcgis <- ggplot(confirmed_cases_linelist) + 
  geom_bar(aes(x=date_confirmation,y=n,fill=pre_reports),stat="identity") + 
  facet_wrap(~province, scales="free_y")
p_confirm_pre_arcgis

####################################
## MERGE ARCGIS AND LINELIST DATA
####################################
## Using linelist data for dates at and before 21.01.2020,
## arcgis data otherwise
use_data_diff$country <- use_data_diff$country_region
combined_dat_final <- merge_data(combined_dat, use_data_diff, switch_date="21.01.2020")
combined_dat_final$ID <- 1:nrow(combined_dat_final)

source("code/fit_delay_distributions.R")



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

####################################
## SYMPTOM ONSET DISTRIBUTION
## UPDATED -- Now using the distribution from Weibull derived by
## Backer et al.
####################################
## 1000 draws from their posterior
n_samps <- 1000
times <- seq(0,25,by=0.1)
weibull_dists <- matrix(0, nrow=1000, ncol=length(times))
for(i in seq_len(n_samps)){
  pars <- weibull_stan_draws[i,]
  alpha <- pars$alpha
  sigma <- pars$sigma
  weibull_dists[i,] <- dweibull(times, alpha, sigma)
}
colnames(weibull_dists) <- times
weibull_dists_bounds <- as.data.frame(t(apply(weibull_dists, 2, function(x) quantile(x, c(0.025,0.5,0.975)))))
colnames(weibull_dists_bounds) <- c("lower","median","upper")
weibull_dists_bounds$times <- times

p_incubation <- ggplot(weibull_dists_bounds) + 
  geom_ribbon(aes(x=times, ymax=upper,ymin=lower),fill="grey70",alpha=0.4,col="black") + 
  geom_line(aes(x=times,y=median),size=1) +
  ylab("Probability density") +
  xlab("Days since onset of infection") +
  ggtitle("Incubation period distribution\n (Weibull, time from infection to symptoms)") +
  scale_y_continuous(limits=c(0,0.3),expand=c(0,0),breaks=seq(0,0.3,by=0.05)) +
  scale_x_continuous(expand=c(0,0)) +
  theme_pubr()
p_incubation


#############################
## FULL AUGMENTATION
#############################

## Now let's repeat this process many times to get a distribution
sim_data_infections <- matrix(NA, nrow=repeats, ncol=nrow(combined_dat_final))
sim_data_symptoms <- matrix(NA, nrow=repeats, ncol=nrow(combined_dat_final))

## For each sample, draw a Weibull distribution from the posterior for 
## the incubation period and generate augmented infection times for all individuals
for(i in seq_len(repeats)){
  ## Random draw from the weibull posterior
  incu_period_rand <- weibull_stan_draws[sample(seq_len(nrow(weibull_stan_draws)),1),]
  alpha <- incu_period_rand$alpha
  sigma <- incu_period_rand$sigma
  # sample from posterior if bayesian
  if(bayesian_p_confirm_delay) {
    p_confirm_delay <- sample(fit_kudos$par,1)
  } else {
    # use point estimate if frequentist
    p_confirm_delay <- fit_kudos$par
  }
  
  ## Get symptom onset and infection times
  tmp <- augment_infection_times(combined_dat_final, 
                                 inc_period_alpha=alpha, 
                                 inc_period_sigma=sigma, 
                                 fit_kudos_gamma$par[1],
                                 fit_kudos_gamma$par[2])
                                 #p_confirm_delay=fit_kudos$par)
  
  sim_data_infections[i,] <- tmp$augmented_infection_times
  sim_data_symptoms[i,] <- tmp$augmented_symptom_onsets
}

sim_data_infections_melted <- reshape2::melt(sim_data_infections)
sim_data_symptoms_melted <- reshape2::melt(sim_data_symptoms)

sim_data_infections_melted$var <- "date_infection"
sim_data_symptoms_melted$var <- "date_onset_symptoms"

colnames(sim_data_infections_melted) <- colnames(sim_data_symptoms_melted) <- c("repeat_no","individual","date","var")

## Combine symptom onsets and infections and conert to dates
sim_data_all <- rbind(sim_data_infections_melted, sim_data_symptoms_melted)
sim_data_all$date <- as.Date(floor(sim_data_all$date), origin="1970-01-01")

## Sum by repeat, variable and date ie. events per day
sim_data_sum <- sim_data_all %>% group_by(repeat_no, var, date) %>% tally()
sim_data_sum <- sim_data_sum %>% group_by(repeat_no, var) %>%
  mutate(date_diff = as.numeric(date_today - date))

## Divide by proportion observed
source("code/unobserved_proportion.R")

sim_data_sum <- sim_data_sum %>% group_by(repeat_no, var) %>%
  mutate(prop_observed = ifelse(var == "date_infection", cumsum(prop_seen)[date_diff], prop_confirmed[date_diff])) #%>%
  #select(-date_diff)
sim_data_sum$n_inflated <- rnbinom(nrow(sim_data_sum), sim_data_sum$n, sim_data_sum$prop_observed) + sim_data_sum$n
#sim_data_sum <- sim_data_sum %>% mutate(n_inflated = ifelse(var=="date_infection",n_inflated, n))
sim_data_sum <- sim_data_sum %>% ungroup() %>% complete(repeat_no, var, date, fill=list(n=0,n_inflated=0,prop_observed=0))

variable_key2 <- c("date_confirmation"="Confirmation date (known)",
                   "date_onset_symptoms"="Number of observed cases with symptom onset (estimated from date of confirmation)",
                   "date_onset_symptoms_inflated"="Number of observed and not yet observed cases with symptom onset (estimated)",
                   "date_admission_hospital"="Hospital admission date",
                   "date_infection"="Number of infections for observed cases (estimated)",
                   "date_infection_inflated" = "Number of infections for observed cases and those not yet observed (estimated) ")

################################################
## OVERALL PLOT
## Distribution of times for each date
sim_data_quantiles <- sim_data_sum %>% group_by(date, var) %>% 
  do(data.frame(t(c(quantile(.$n, probs = c(0.025,0.5,0.975),na.rm=TRUE),mean(.$n)))))

sim_data_quantiles_inflated <- sim_data_sum %>% group_by(date, var) %>% 
  do(data.frame(t(c(quantile(.$n_inflated, probs = c(0.025,0.5,0.975),na.rm=TRUE),mean(.$n_inflated)))))

sim_data_quantiles_inflated$var <- c("date_infection" = "date_infection_inflated",
                                     "date_onset_symptoms" = "date_onset_symptoms_inflated")[sim_data_quantiles$var]

sim_data_quantiles <- bind_rows(sim_data_quantiles, sim_data_quantiles_inflated) %>% arrange(date)
## Get confirmation time data
confirm_data <- combined_dat_final %>% filter(!is.na(date_confirmation)) %>% group_by(date_confirmation) %>% tally()
confirm_data$Variable <- "Number of confirmations for observed cases (observed)"

sim_data_quantiles$var <- variable_key2[sim_data_quantiles$var]

colnames(sim_data_quantiles) <- c("date","Variable","lower","median","upper","mean")

tmp <- which(rev(cumsum(prop_seen)) > 0.95)
tmp[length(tmp)]
threshold_99 <- convert_date(date_today) + times[tmp[length(tmp)]]

tmp <- which(rev(cumsum(prop_seen)) > 0.8)
tmp[length(tmp)]
threshold_80 <- convert_date(date_today) + times[tmp[length(tmp)]]

tmp <- which(rev(cumsum(prop_seen)) > 0.5)
tmp[length(tmp)]
threshold_50 <- convert_date(date_today) + times[tmp[length(tmp)]]

tmp <- which(rev(cumsum(prop_seen)) > 0.2)
tmp[length(tmp)]
threshold_20 <- convert_date(date_today) + times[tmp[length(tmp)]]

thresholds <- c(threshold_99, threshold_80, threshold_50, threshold_20)


#augmented_data_plot <- plot_augmented_data(sim_data_quantiles, confirm_data,ymax=2000,ybreaks=100,max_date = date_today, thresholds)

sim_data_quantiles_truncated <- sim_data_quantiles %>% filter(date <= convert_date(Sys.Date()))
augmented_data_plot <- plot_augmented_data(sim_data_quantiles_truncated, confirm_data,ymax=6000,ybreaks=500,
                                           max_date = date_today, min_date="01.01.2020", thresholds=thresholds)
augmented_data_plot

onset_only <- sim_data_quantiles_truncated %>% 
                filter(Variable %in% c("Number of observed cases with symptom onset (estimated from date of confirmation)",
                                  "Number of observed and not yet observed cases with symptom onset (estimated)"))
augmented_plot_onset <- plot_augmented_data(onset_only, confirm_data,ymax=5000,ybreaks=500,
                                            max_date = date_today, min_date="01.01.2020", thresholds=NULL,
                                            cols = c("grey40","orange","red"), cols2 = c("orange","red"),
                                            title = "Augmented and observed timings of symptom onset in China")
augmented_plot_onset
infection_only <- sim_data_quantiles_truncated %>% 
  filter(Variable %in% c("Number of infections for observed cases (estimated)",
                         "Number of infections for observed cases and those not yet observed (estimated) "))
augmented_plot_infection <- plot_augmented_data(infection_only, confirm_data,ymax=5000,ybreaks=500,
                                                max_date = date_today, min_date="01.01.2020", thresholds=NULL,
                                                cols = c("grey40","blue","skyblue"), cols2 = c("blue","skyblue"),
                                                title = "Augmented and observed timings of infection in China")
augmented_plot_infection
#plot_grid(augmented_plot_onset, augmented_plot_infection)
## Distribution of times for each individual

sim_data_quantiles_indiv <- sim_data_all %>% group_by(individual, var) %>% 
  do(data.frame(t(quantile(as.numeric(.$date), probs = c(0.025,0.5,0.975),na.rm=TRUE))))

sim_data_quantiles_indiv$X2.5. <- as.Date(sim_data_quantiles_indiv$X2.5., origin="1970-01-01")
sim_data_quantiles_indiv$X50. <- as.Date(sim_data_quantiles_indiv$X50., origin="1970-01-01")
sim_data_quantiles_indiv$X97.5. <- as.Date(sim_data_quantiles_indiv$X97.5., origin="1970-01-01")


#########################
## FINAL HOUSEKEEPING
## Tidy up data to share
sim_data_infections1 <- as.data.frame(t(sim_data_infections[1:100,]))
for(i in seq_len(ncol(sim_data_infections1))){
  sim_data_infections1[,i] <- as.Date(floor(sim_data_infections1[,i]), origin="1970-01-01")
}
sim_data_infections1 <- bind_cols(combined_dat_final,sim_data_infections1)


sim_data_symptoms1 <- as.data.frame(t(sim_data_symptoms[1:100,]))
for(i in seq_len(ncol(sim_data_symptoms1))){
  sim_data_symptoms1[,i] <- as.Date(floor(sim_data_symptoms1[,i]), origin="1970-01-01")
}
sim_data_infections1 <- bind_cols(combined_dat_final,sim_data_symptoms1)

write_csv(sim_data_infections1, path="augmented_data/augmented_infection_times.csv")
write_csv(sim_data_symptoms1, path="augmented_data/augmented_symptom_times.csv")

## Create results panel plot programmatically
element_text_size <- 11
text_size_theme <- theme(title=element_text(size=element_text_size), 
                         axis.text=element_text(size=element_text_size), 
                         axis.title = element_text(size=element_text_size))
p_other_confirm_fit1 <- p_other_confirm_fit + text_size_theme
p_incubation1 <- p_incubation + text_size_theme
p_confirm_delay_kudos <- p_confirm_delay_kudos + text_size_theme
assumption_plot <- plot_grid(p_confirm_delay_kudos, p_incubation1,ncol=2,align="hv")
augmented_data_plot1 <- augmented_data_plot + theme(legend.position=c(0.2,0.2)) 
layout <- c(
  area(t=0,b=12,l=0,r=18),
  area(t=2,b=7,l=2,r=14)
)

results_panel <- augmented_data_plot1 + assumption_plot + plot_layout(design=layout)
results_panel
#######################
## SPATIAL PLOTS
#######################
individual_key <- combined_dat_final[,c("ID","age","country","sex","city","province","latitude","longitude")]
colnames(individual_key)[1] <- "individual"
sim_data_all <- as_tibble(sim_data_all)

merged_data <- right_join(individual_key, sim_data_all,by=c("individual"))
merged_data <- merged_data[!is.na(merged_data$date),]

start_dates <- merged_data %>% group_by(var, province,repeat_no) %>% 
  mutate(first_day=min(date)) %>% ungroup() %>% 
  group_by(var, province) %>%
  filter(var=="date_infection") 

ggplot(start_dates) +
  geom_histogram(aes(x=first_day)) + facet_wrap(~province,scales="free_y")

#############################
## Aggregate by province
## Get confirmation time data
confirm_dat_province <- combined_dat_final %>% ungroup() %>% filter(!is.na(combined_dat_final$date_confirmation)) %>% 
  group_by(province, date_confirmation) %>% tally() %>%
  ungroup() %>% complete(province, date_confirmation, fill=list(n=0))
confirm_dat_province$Variable <- "Confirmed cases of infections that have been observed"

## for each province, variable and sample,
## Find the first date of an infection
## Then, get the average across all repeats
## Shift all dates so that this date is day 0
#merged_data <- merged_data %>% 
#  group_by(province, var, repeat_no) %>% 
#  mutate(date = ifelse(date < convert_date("01.11.2019"),convert_date("01.11.2019"),date)) %>%
#  mutate(start_day=min(date)) %>% ungroup() %>%
#  group_by(province, var,repeat_no) %>%
#  mutate(d_diff_mean=as.numeric(date - start_day))



province_data <- merged_data %>% 
  group_by(repeat_no, var, date, province) %>%
  tally() 

province_data <- province_data %>% group_by(repeat_no, var, province) %>%
  mutate(date_diff = as.numeric(date_today - date))
province_data <- province_data %>% group_by(repeat_no, var, province) %>%
  mutate(prop_observed = ifelse(var == "date_infection", cumsum(prop_seen)[date_diff], prop_confirmed[date_diff]),
         n_inflated=floor(n/prop_observed))# %>%
#select(-date_diff)
province_data$n_inflated <- rnbinom(nrow(province_data), province_data$n, province_data$prop_observed) + province_data$n
#sim_data_sum <- sim_data_sum %>% mutate(n_inflated = ifelse(var=="date_infection",n_inflated, n))
province_data <- province_data %>% ungroup() %>% complete(repeat_no, var, date, province, fill=list(n=0,n_inflated=0,prop_observed=0))


sim_data_quantiles_province <- province_data %>% group_by(date, var, province) %>% 
  do(data.frame(t(quantile(.$n_inflated, probs = c(0.025,0.5,0.975),na.rm=TRUE))))

sim_data_quantiles_province$var <- variable_key2[sim_data_quantiles_province$var]
colnames(sim_data_quantiles_province) <- c("date","Variable","province","lower","median","upper")
total_confirmed_prov <- confirm_dat_province %>% group_by(province) %>% summarise(n=sum(n))
total_confirmed_prov <- total_confirmed_prov[order(-total_confirmed_prov$n),]
factor_order <- as.character(total_confirmed_prov$province)

confirm_dat_province$province <- factor(as.character(confirm_dat_province$province), 
                                         levels=factor_order)
sim_data_quantiles_province$province <- factor(as.character(sim_data_quantiles_province$province), 
                                               levels=factor_order)


by_province <- plot_augmented_data_province(sim_data_quantiles_province, confirm_dat_province)
top_6 <- factor_order[1:6]
by_province_top6 <- plot_augmented_data_province(sim_data_quantiles_province[sim_data_quantiles_province$province %in% top_6,], 
                                                 confirm_dat_province[confirm_dat_province$province %in% top_6,], max_date=date_today)
by_province_top6 <- by_province_top6 + facet_wrap(~province, ncol=3, scales="free_y") + 
  theme(legend.text=element_text(size=10))
by_province_top6

## Plot time from start
p_start_delay_dist <- plot_time_from_start(sim_data_infections_melted, individual_key,xmax=100)
p_start_delay_dist

# source("code/shifting_curves.R")
