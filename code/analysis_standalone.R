#setwd("~/Documents/case_to_infection/")
setwd("~/GitHub/case_to_infection/")
refit_p_confirm_delay <- FALSE # if TRUE, fit geometric distribution to confirmation delay data;
# if FALSE, read from file
bayesian_p_confirm_delay <- FALSE # if TRUE, use posterior for confirmation delay parameter, if FALSE, use point estimate

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
source("code/date_functions.R")
source("code/augmentation_functions.R")


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


#############################
## FULL AUGMENTATION
#############################

## Now let's repeat this process many times to get a distribution
sim_data_infections <- matrix(NA, nrow=repeats, ncol=nrow(combined_dat_final))
sim_data_symptoms <- matrix(NA, nrow=repeats, ncol=nrow(combined_dat_final))
augmentation_tracker <- NULL

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
                                 gamma_mean=fit_kudos_gamma$par[1],
                                 gamma_var=fit_kudos_gamma$par[2],
                                 p_confirm_delay=fit_kudos$par)
  
  sim_data_infections[i,] <- tmp$augmented_infection_times
  sim_data_symptoms[i,] <- tmp$augmented_symptom_onsets
  augmentation_tracker[[i]] <- c(i,incu_period_rand$alpha, incu_period_rand$sigma)
  
}
used_weibull_pars <- as_tibble(do.call("rbind", augmentation_tracker))
colnames(used_weibull_pars) <- c("repeat_no","alpha","sigma")


sim_data_infections_melted <- as_tibble(reshape2::melt(sim_data_infections))
sim_data_symptoms_melted <- as_tibble(reshape2::melt(sim_data_symptoms))

sim_data_infections_melted$var <- "date_infection"
sim_data_symptoms_melted$var <- "date_onset_symptoms"

colnames(sim_data_infections_melted) <- colnames(sim_data_symptoms_melted) <- c("repeat_no","individual","date","var")

## Combine symptom onsets and infections and conert to dates
sim_data_all <- rbind(sim_data_infections_melted, sim_data_symptoms_melted)
sim_data_all$date <- as.Date(floor(sim_data_all$date), origin="1970-01-01")
sim_data_all <- as_tibble(sim_data_all)
sim_data_all <- sim_data_all %>% pivot_wider(values_from=date,names_from=var)
sim_data_all <- left_join(sim_data_all, used_weibull_pars)

sim_data_all <- sim_data_all %>% mutate(symp_delay=as.numeric(date_today-date_onset_symptoms),
                                        confirm_delay=as.numeric(date_today-date_infection),
                                        total_delay=symp_delay+confirm_delay)

all_probs_forward <- generate_forward_probabilities_dist(repeats, used_weibull_pars, fit_kudos$par,tmax=ceiling(max(sim_data_all$total_delay)))

sim_data_all <- sim_data_all %>% left_join(all_probs_forward)
sim_data_all <- sim_data_all %>% mutate(n_inflated_inf = rnbinom(n(),1,cumu_prob_total)+1,
                                        n_inflated_symp = rnbinom(n(),1,cumu_prob_confirm)+1)

## Sum by repeat, variable and date ie. events per day
sim_inf_sum <- sim_data_all %>% group_by(repeat_no, date_infection) %>% tally()
colnames(sim_inf_sum)[2] <- "date"
sim_inf_sum$var <- "date_infection"
sim_symp_sum <- sim_data_all %>% group_by(repeat_no, date_onset_symptoms) %>% tally()
colnames(sim_symp_sum)[2] <- "date"
sim_symp_sum$var <- "date_onset_symptoms"

sim_data_sum <- bind_rows(sim_inf_sum, sim_symp_sum)

sim_inf_sum_inflated <- sim_data_all %>% group_by(repeat_no, date_infection) %>% summarise(n_inflated=sum(n_inflated_inf))
colnames(sim_inf_sum_inflated)[2] <- "date"
sim_inf_sum_inflated$var <- "date_infection"
sim_symp_sum_inflated <- sim_data_all %>% group_by(repeat_no, date_onset_symptoms) %>% summarise(n_inflated=sum(n_inflated_symp))
colnames(sim_symp_sum_inflated)[2] <- "date"
sim_symp_sum_inflated$var <- "date_onset_symptoms"

sim_data_sum_inflated <- bind_rows(sim_inf_sum_inflated, sim_symp_sum_inflated)

sim_data_sum <- left_join(sim_data_sum, sim_data_sum_inflated)

sim_data_sum1 <- sim_data_sum %>% complete(repeat_no, var, date, fill=list(n=0,n_inflated=0))

variable_key2 <- c("date_confirmation"="Confirmation date (known)",
                   "date_onset_symptoms"="Number of observed cases with symptom onset (estimated from date of confirmation)",
                   "date_onset_symptoms_inflated"="Number of observed and not yet observed cases with symptom onset (estimated)",
                   "date_admission_hospital"="Hospital admission date",
                   "date_infection"="Number of infections for observed cases (estimated)",
                   "date_infection_inflated" = "Number of infections for observed cases and those not yet observed (estimated)")

################################################
## OVERALL PLOT
## Distribution of times for each date
sim_data_sum_all <- sim_data_sum1 %>% pivot_longer(cols=c("n", "n_inflated"), names_to="inflated",values_to="value")

sim_data_quantiles <- sim_data_sum_all %>% group_by(date, var, inflated) %>% 
  do(data.frame(t(c(quantile(.$value, probs = c(0.025,0.5,0.975),na.rm=TRUE),mean(.$value)))))

sim_data_quantiles_inflated <- sim_data_sum %>% group_by(date, var) %>% 
  do(data.frame(t(c(quantile(.$n_inflated, probs = c(0.025,0.5,0.975),na.rm=TRUE),mean(.$n_inflated)))))

sim_data_quantiles_inflated$var <- c("date_infection" = "date_infection_inflated",
                                     "date_onset_symptoms" = "date_onset_symptoms_inflated")[sim_data_quantiles$var]

sim_data_quantiles <- bind_rows(sim_data_quantiles, sim_data_quantiles_inflated) %>% arrange(date)
## Get confirmation time data
confirm_data <- combined_dat_final %>% filter(!is.na(date_confirmation)) %>% group_by(date_confirmation) %>% tally()
confirm_data$Variable <- "Number of confirmations for observed cases (observed)"

sim_data_quantiles$var <- variable_key2[sim_data_quantiles$var]

colnames(sim_data_quantiles) <- c("date","Variable","inflated","lower","median","upper","mean")

tmp <- which(rev(cumsum(prop_seen)) > 0.95)
#tmp[length(tmp)]
threshold_95 <- convert_date(date_today) + times[tmp[length(tmp)]]

tmp <- which(rev(cumsum(prop_seen)) > 0.8)
#tmp[length(tmp)]
threshold_80 <- convert_date(date_today) + times[tmp[length(tmp)]]

tmp <- which(rev(cumsum(prop_seen)) > 0.5)
#tmp[length(tmp)]
threshold_50 <- convert_date(date_today) + times[tmp[length(tmp)]]

tmp <- which(rev(cumsum(prop_seen)) > 0.2)
#tmp[length(tmp)]
threshold_20 <- convert_date(date_today) + times[tmp[length(tmp)]]

thresholds <- c(threshold_95, threshold_80, threshold_50, threshold_20)


#augmented_data_plot <- plot_augmented_data(sim_data_quantiles, confirm_data,ymax=2000,ybreaks=100,max_date = date_today, thresholds)

sim_data_quantiles_truncated <- sim_data_quantiles %>% filter(date <= convert_date(Sys.Date()))
augmented_data_plot <- plot_augmented_data(sim_data_quantiles_truncated, confirm_data,ymax=6000,ybreaks=500,
                                           max_date = date_today, min_date="01.01.2020", thresholds=NULL)
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
                         "Number of infections for observed cases and those not yet observed (estimated)"))
infection_only[infection_only$Variable == "Number of infections for observed cases and those not yet observed (estimated)","mean"] <- NA

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
