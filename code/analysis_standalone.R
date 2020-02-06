######################
## SETUP
setwd("~/Documents/case_to_infection/")
#setwd("~/GitHub/case_to_infection/")
savewd <- "plots1"
refit_p_confirm_delay <- TRUE # if TRUE, fit geometric distribution to confirmation delay data;
# if FALSE, read from file
bayesian_p_confirm_delay <- FALSE # if TRUE, use posterior for confirmation delay parameter, if FALSE, use point estimate

use_geometric_confirmation_delay <- FALSE
save_augmented_results <- FALSE

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
library(multidplyr)

## For use with multidplyr
cluster <- new_cluster(4)
cluster_library(cluster, "tidyverse")

## Some setup for STAN
if(refit_p_confirm_delay && bayesian_p_confirm_delay) {
  library(rstan)
  options(mc.cores = parallel::detectCores())
  rstan_options(auto_write = TRUE)
}
source("code/analysis_functions.R")
source("code/date_functions.R")
source("code/augmentation_functions.R")

## Need to be careful here - today's date needs to be
## the last day at which there are final case counts
## for that day
date_today <- convert_date("04.02.2020")

weibull_stan_draws <- read.csv("data/backer_weibull_draws.csv")
minimum_confirmation_delay <- 1
## source("code/plot_china_map.R")

## These column names will be kept as keys for the rest of the analysis
key_colnames <- c("ID","age","country", "sex","city","province",
                  "latitude","longitude")
## These are the column names with used variables
var_colnames <- c("date_confirmation","date_onset_symptoms","date_admission_hospital")
use_colnames <- c(key_colnames, var_colnames)

## Number of bootstrap samples to take. Set this to something small for a quick run
repeats <- 100

#########################
## LOAD DATA
## load the data - try to only do this once otherwise auth token gets stale
## First step is to clean and take a look at the data
## This combines the data for Hubei and other locations in China
## Note this is ONLY China

kudos_dat <- readRDS("data/kudos_dat.rds")
combined_dat <- readRDS("data/combined_dat.rds")
use_data_diff <- readRDS("data/use_data_diff.rds")
#use_data_subset <- readRDS("data/use_data_subset.rds")

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

kudos_dat_china <- kudos_dat %>% mutate(hosp_delay = hosp_visit_date - symptom_onset) %>% 
  filter(country %in% c("China","Taiwan","Hong Kong"))

 saveRDS(kudos_dat, "data/kudos_dat.rds")
 saveRDS(combined_dat, "data/combined_dat.rds")
 #saveRDS(use_data_subset, "data/use_data_subset.rds")
 saveRDS(use_data_diff, "data/use_data_diff.rds")

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

## Only using line list data up to and including 21.01.2020
confirmed_cases_linelist <- confirmed_cases_linelist %>% 
  mutate(pre_reports=date_confirmation <= convert_date("21.01.2020"))

## Have a quick look to see how this compares to the ARCGIS data
p_confirm_pre_arcgis <- ggplot(confirmed_cases_linelist) + 
  geom_bar(aes(x=date_confirmation,y=n,fill=pre_reports),stat="identity") + 
  facet_wrap(~province, scales="free_y")

####################################
## MERGE ARCGIS AND LINELIST DATA
####################################
## Using linelist data for dates at and before 21.01.2020,
## arcgis data otherwise

## Death data
combined_dat_death <- combined_dat %>% filter(outcome == "died" & !is.na(date_death_or_discharge))


use_data_diff$country <- use_data_diff$country_region
combined_dat_final <- merge_data(combined_dat, use_data_diff, switch_date="21.01.2020")
combined_dat_final <- combined_dat_final %>% filter(date_confirmation <= date_today)
combined_dat_final$individual <- 1:nrow(combined_dat_final)

source("code/fit_delay_distributions.R")
source("code/generate_delay_distributions.R")
## Alternative using the  other line list data
## source("code/generate_delay_distributions_moritz_data.R")
## Create results panel plot programmatically
element_text_size <- 11
text_size_theme <- theme(title=element_text(size=element_text_size), 
                         axis.text=element_text(size=element_text_size), 
                         axis.title = element_text(size=element_text_size))
p_other_confirm_fit1 <- p_other_confirm_fit + text_size_theme
p_incubation1 <- p_incubation + text_size_theme
p_confirm_delay_kudos <- p_confirm_delay_kudos + text_size_theme
assumption_plot <- p_incubation1 / p_confirm_delay_kudos
assumption_plot

p_sliding_delays_forward
p_sliding_delays_backward

## Fill up backward confirmation delay distribution to first confirmation date
min_date <- min(combined_dat_final$date_confirmation)
dates <- convert_date(min_date:(min(gamma_pars_dat_backward$date_confirmation)-1))
gamma_mean_use <- gamma_pars_dat_backward %>% filter(date_confirmation == min(date_confirmation)) %>% pull(gamma_mean_backward)
gamma_var_use <- gamma_pars_dat_backward %>% filter(date_confirmation == min(date_confirmation)) %>% pull(gamma_var_backward)
gamma_pars_dat_backward <- gamma_pars_dat_backward %>% bind_rows(tibble(date_confirmation=dates, gamma_mean_backward=gamma_mean_use, 
                                                       gamma_var_backward=gamma_var_use, n_used=threshold, direction="backward")) %>% arrange(date_confirmation) %>% as_tibble
gamma_pars_dat_backward <- gamma_pars_dat_backward  %>% mutate(gamma_scale_backward=gamma_var_backward/gamma_mean_backward,
                                             gamma_shape_backward=gamma_mean_backward/gamma_scale_backward) %>%
  select(-c("gamma_mean_backward","gamma_var_backward"))
combined_dat_final <- combined_dat_final %>% left_join(gamma_pars_dat_backward)

#############################
## FULL AUGMENTATION
#############################
## Now let's repeat this process many times to get a distribution
sim_data_infections <- matrix(NA, ncol=repeats, nrow=nrow(combined_dat_final))
sim_data_symptoms <- matrix(NA, ncol=repeats, nrow=nrow(combined_dat_final))
augmentation_tracker <- NULL

## For each sample, draw a Weibull distribution from the posterior for 
## the incubation period and generate augmented infection times for all individuals
for(i in seq_len(repeats)){
  ## Random draw from the weibull posterior
  incu_period_rand <- weibull_stan_draws[sample(seq_len(nrow(weibull_stan_draws)),1),]
  alpha <- incu_period_rand$alpha
  sigma <- incu_period_rand$sigma
  
  p_confirm_delay <- NULL
  # sample from posterior if bayesian
  if(use_geometric_confirmation_delay) {
    if(bayesian_p_confirm_delay) {
      p_confirm_delay <- sample(fit_kudos$par,1)
    } else {
      # use point estimate if frequentist
      p_confirm_delay <- fit_kudos$par
    }
  } 
    
  ## Get symptom onset and infection times
  tmp <- augment_infection_times(combined_dat_final, 
                                 inc_period_alpha=alpha, 
                                 inc_period_sigma=sigma, 
                                 p_confirm_delay=p_confirm_delay,
                                 minimum_confirmation_delay)
  data_date_confirmation <- combined_dat_final$date_confirmation
  individuals <- combined_dat_final$individual
  sim_data_infections[,i] <- tmp$augmented_infection_times
  sim_data_symptoms[,i] <- tmp$augmented_symptom_onsets
  
  ## Keep track of which incubation period parameters we used
  augmentation_tracker[[i]] <- c(i,incu_period_rand$alpha, incu_period_rand$sigma)
}
#######################################
## AUGMENTATION DONE, NOW HOUSEKEEPING
#######################################
## Put outputs into nice tibbles
used_weibull_pars <- as_tibble(do.call("rbind", augmentation_tracker))
colnames(used_weibull_pars) <- c("repeat_no","alpha","sigma")

sim_data_infections <- data.frame(sim_data_infections)
sim_data_infections$individual <- individuals
sim_data_infections$date_confirmation <- data_date_confirmation

sim_data_symptoms <- data.frame(sim_data_symptoms)
sim_data_symptoms$individual <- individuals
sim_data_symptoms$date_confirmation <- data_date_confirmation

sim_data_infections_melted <- as_tibble(sim_data_infections) %>% 
  pivot_longer(-c("date_confirmation", "individual"), names_to="repeat_no",values_to="date_infection") %>%
  mutate(repeat_no=match(repeat_no, repeat_no))
sim_data_symptoms_melted <-as_tibble(sim_data_symptoms) %>% 
  pivot_longer(-c("date_confirmation", "individual"), names_to="repeat_no",values_to="date_onset_symptoms") %>%
  mutate(repeat_no=match(repeat_no, repeat_no))

sim_data_all <- left_join(sim_data_infections_melted, sim_data_symptoms_melted)

## Combine symptom onsets and infections and convert to dates
sim_data_all$date_infection <- as.Date(floor(sim_data_all$date_infection), origin="1970-01-01")
sim_data_all$date_onset_symptoms <- as.Date(floor(sim_data_all$date_onset_symptoms), origin="1970-01-01")

## Fill up FORWARD confirmation delay distribution to first confirmation date
min_date <- min(sim_data_all$date_onset_symptoms)
dates <- convert_date(min_date:(min(gamma_pars_dat_forward$date_onset_symptoms)-1))
gamma_mean_use <- gamma_pars_dat_forward %>% filter(date_onset_symptoms == min(date_onset_symptoms)) %>% pull(gamma_mean_forward)
gamma_var_use <- gamma_pars_dat_forward %>% filter(date_onset_symptoms == min(date_onset_symptoms)) %>% pull(gamma_var_forward)
gamma_pars_dat_forward <- gamma_pars_dat_forward %>% bind_rows(tibble(date_onset_symptoms=dates, gamma_mean_forward=gamma_mean_use, 
                                                                        gamma_var_forward=gamma_var_use, n_used=threshold,direction="forward")) %>% 
  arrange(date_onset_symptoms) %>% as_tibble
gamma_pars_dat_forward <- gamma_pars_dat_forward  %>% mutate(gamma_scale_forward=gamma_var_forward/gamma_mean_forward,
                                                               gamma_shape_forward=gamma_mean_forward/gamma_scale_forward) %>%
  select(-c("gamma_mean_forward","gamma_var_forward","direction"))


## Get parameters for nowcasting stuff
sim_data_all <- left_join(sim_data_all, used_weibull_pars) %>% left_join(gamma_pars_dat_forward)
sim_data_all <- sim_data_all %>% mutate(fixed_geom=fit_kudos$par)

## Get delays for all augmented and actual events
## i) delay from infection to symptom onset
## ii delay from symptom onset to confirmation
## iii) total delay
sim_data_all <- sim_data_all %>% mutate(symp_delay=as.numeric(date_onset_symptoms-date_infection),
                                        confirm_delay=as.numeric(date_confirmation-date_onset_symptoms),
                                        total_delay=symp_delay+confirm_delay)

## Free some memory
rm(sim_data_symptoms_melted)
rm(sim_data_infections_melted)
rm(sim_data_infections)
rm(sim_data_symptoms)

all_delay_prob_parameters <- sim_data_all %>% 
  select(repeat_no, alpha, sigma, gamma_scale_forward, gamma_shape_forward, fixed_geom, date_confirmation, date_onset_symptoms) %>% 
  distinct()

## Get table of probabilities that events have happened after certain delays
## This takes a while
if(!exists("all_probs_forward")) {
  all_probs_forward <- generate_forward_probabilities_dist(repeats, all_delay_prob_parameters, tmax=300)
}
confirm_probs_gamma <- all_probs_forward[[1]] %>% select(repeat_no, date_onset_symptoms, confirm_delay, cumu_prob_confirm)
colnames(confirm_probs_gamma)[2] <- "date"
confirm_probs_geometric <- all_probs_forward[[2]] %>% select(repeat_no, confirm_delay, cumu_prob_confirm)
symptom_probs <- all_probs_forward[[3]] %>% select(repeat_no, symp_delay, cumu_prob_symp)

source("code/generate_inflations.R")

rm(symptom_all)
rm(infections_all)

################################################
## OVERALL PLOT
## Distribution of times for each date
final_quantiles <- final_all %>% select(repeat_no, var, date, inflated0, total) %>% 
  pivot_longer(cols=c("inflated0","total"),names_to="inflated") %>% 
  group_by(date, var, inflated) %>% 
  do(data.frame(t(c(quantile(.$value, probs = c(0.01,0.025,0.25,0.5,0.75,0.975,0.99),na.rm=TRUE),mean(.$value)))))
colnames(final_quantiles) <- c("date","var","inflated",
                               "min","lower","midlow","median","midhigh","upper","max","mean")
## Get confirmation time data
confirm_data <- combined_dat_final %>% filter(!is.na(date_confirmation)) %>% group_by(date_confirmation) %>% tally()
confirm_data$inflated <- "total"
confirm_data$var <- "confirmed"

## Get the same but cumulative
final_all_cumulative <- final_all %>% group_by(repeat_no, var) %>% arrange(repeat_no, var, date) %>% 
  mutate(cumu_observed=cumsum(inflated0),
         cumu_unobserved=cumsum(inflated1),
         cumu_total=cumsum(total)) %>% ungroup()
final_quantiles_cumulative <- final_all_cumulative %>% select(repeat_no, var, date, cumu_observed, cumu_total) %>% 
  pivot_longer(cols=c("cumu_observed","cumu_total"),names_to="inflated") %>% 
  group_by(date, var, inflated) %>% 
  do(data.frame(t(c(quantile(.$value, probs = c(0.01,0.025,0.25,0.5,0.75,0.975,0.99),na.rm=TRUE),mean(.$value)))))
colnames(final_quantiles_cumulative) <- c("date","var","inflated",
                               "min","lower","midlow","median","midhigh","upper","max","mean")
change_key <- c("cumu_observed"="inflated0", "cumu_total"="total")
final_quantiles_cumulative$inflated <- change_key[final_quantiles_cumulative$inflated]

## Get confirmation time data
confirm_data_cumulative <- combined_dat_final %>% filter(!is.na(date_confirmation)) %>% 
  group_by(date_confirmation) %>% tally() %>% mutate(n=cumsum(n)) %>% ungroup()
confirm_data_cumulative$inflated <- "total"
confirm_data_cumulative$var <- "confirmed"


## Get vertical dashes to show confirmation proportions over time
prop_seen_mean <- generate_total_forward_probabilities_dist(repeats=repeats, all_delay_prob_parameters,
                                                       tmax=100,use_geometric_confirmation_delay1=use_geometric_confirmation_delay)
times <- rev(seq(date_today-200, date_today,by=1))

rm(all_probs_forward)

###################################################
## GET PROPORTION OBSERVED BY DATE OF INFECTION
## Get times at which at least these % of infections from that day
## have been seen
threshold_vals <- c(0.95, 0.8, 0.5, 0.2)
if (use_geometric_confirmation_delay) {
  prop_sympt_observed_mean <- confirm_probs_geometric %>% group_by(confirm_delay) %>% summarise(cumu_prob_total=mean(cumu_prob_confirm))
  ## Get mean overall confirmation probs
  prop_seen <- prop_seen_mean %>% pull(cumu_prob_total)
  ## Sequence from today to 200 days ago
  ## First day in the past at which these percentages have been seen
  thresholds <- times[sapply(threshold_vals, function(x) which(prop_seen > x)[1])]
} else {
  prop_sympt_observed_mean <- confirm_probs_gamma %>% 
    mutate(actual_delay=as.numeric(date_today-date-1)) %>% 
    filter(actual_delay == confirm_delay) %>% group_by(confirm_delay, date) %>% 
    summarise(cumu_prob_total=mean(cumu_prob_confirm))
  
  
  thresholds <- prop_seen_mean %>% mutate(seen_95=ifelse(cumu_prob_total > 0.95,1,0),
                            seen_80=ifelse(cumu_prob_total > 0.80,1,0),
                            seen_50=ifelse(cumu_prob_total > 0.50,1,0),
                            seen_20=ifelse(cumu_prob_total > 0.20,1,0)) %>% 
    pivot_longer(c("seen_95","seen_80","seen_50","seen_20"),
                   names_to="threshold",values_to="passed") %>%
    filter(passed==1) %>% group_by(threshold) %>% summarise(res=max(date_infection)) %>% pull(res) %>% rev()
  
}
#####################################################
#####################################################
## GET PROPORTION OBSERVED BY DATE OF SYMPTOM ONSET
prop_symp_seen <- prop_sympt_observed_mean %>% pull(cumu_prob_total)
thresholds_symp <- times[sapply(threshold_vals, function(x) which(prop_symp_seen > x)[1])]

p_result <- plot_augmented_data(final_quantiles, confirm_data, max_date=date_today, min_date="15.12.2019",
                          ymax1=5000,ymax2=5000,ybreaks=1000,thresholds=thresholds,thresholds_symp = thresholds_symp)
p_result_cumu <- plot_augmented_data(final_quantiles_cumulative, confirm_data_cumulative, max_date=date_today, min_date="15.12.2019",
                                ymax1=50000,ymax2=30000,ybreaks=5000,thresholds=thresholds,thresholds_symp = thresholds_symp)


rm(symptom_observed)
rm(symptom_unobserved)
rm(infections_unobserved)
rm(final_infections_tally)
rm(final_symptom_onsets_tally)
gc()

#######################
## SPATIAL PLOTS
#######################
#############################
## Aggregate by province
## Get confirmation time data
## Shoved it in a script because it's long...
source("code/generate_byprovince_inflations.R")

rm(sim_data_all)

rm(symptom_all_province)
rm(infections_all_province)


final_quantiles_province <- final_all_province %>% select(repeat_no, var, date, inflated0, total, province) %>% 
  pivot_longer(cols=c("inflated0","total"),names_to="inflated") %>% 
  group_by(date, var, inflated, province) %>% 
  do(data.frame(t(c(quantile(.$value, probs = c(0.01,0.025,0.25,0.5,0.75,0.975,0.99),na.rm=TRUE),mean(.$value)))))


colnames(final_quantiles_province) <- c("date","var","inflated","province",
                                        "min","lower","midlow","median",
                                        "midhigh","upper","max","mean")
final_quantiles$var_full <- paste0(final_quantiles$var, "_", final_quantiles$inflated)

confirm_dat_province <- combined_dat_final %>% ungroup() %>% filter(!is.na(combined_dat_final$date_confirmation)) %>% 
  group_by(province, date_confirmation) %>% tally() %>%
  ungroup() %>% complete(province, date_confirmation, fill=list(n=0))
confirm_dat_province$Variable <- "Confirmed cases of infections that have been observed"
confirm_dat_province$inflated <- "total"
confirm_dat_province$var <- "confirmed"


total_confirmed_prov <- confirm_dat_province %>% group_by(province) %>% summarise(n=sum(n))
total_confirmed_prov <- total_confirmed_prov[order(-total_confirmed_prov$n),]
factor_order <- as.character(total_confirmed_prov$province)

confirm_dat_province$province <- factor(as.character(confirm_dat_province$province), 
                                         levels=factor_order)
final_quantiles_province$province <- factor(as.character(final_quantiles_province$province), 
                                               levels=factor_order)
top_6_provinces <- factor_order[1:6]

p_infections <- plot_augmented_events_byprovince(data_quantiles_province=final_quantiles_province, 
                                 confirmed_data_province=confirm_dat_province,
                                 provinces = top_6_provinces,
                                 var_name="date_infections",max_date="04.02.2020",min_date="01.12.2019",
                                 thresholds=NULL,ncol=3)
p_symptoms <- plot_augmented_events_byprovince(data_quantiles_province=final_quantiles_province, 
                                                 confirmed_data_province=confirm_dat_province,
                                               provinces=top_6_provinces,
                                                 var_name="date_onset_symptoms",max_date="03.02.2020",min_date="01.12.2019",
                                                 cols=c("orange","red"),
                                                 thresholds=NULL, ncol=3)

## Plot time from start
#p_start_delay_dist <- plot_time_from_start(sim_data_infections_melted, individual_key,xmax=100)
#p_start_delay_dist

source("code/shifting_curves.R")


#########################
## FINAL HOUSEKEEPING
## Tidy up data to share
if (save_augmented_results) {
  n_subset <- 100
  all_repeats <- 1:repeats
  use_repeats <- sample(all_repeats, n_subset)
  final_dat_to_share <- final_all %>% filter(repeat_no %in% use_repeats)
  colnames(final_dat_to_share) <- c("repeat_no","variable","date","from_confirmed","inflated_event","total_events")
  
  final_infections_share <- infections_all %>% select(repeat_no, individual, date_infection, symp_delay, augmented) %>% 
    filter(repeat_no %in% use_repeats)
  final_symptom_onsets_share <- symptom_all %>% select(repeat_no, individual, date_infection, 
                                                       date_onset_symptoms, symp_delay, confirm_delay, total_delay, augmented) %>%
    filter(repeat_no %in% use_repeats)
  
  write_csv(final_dat_to_share, path="augmented_data/augmented_totals.csv")
  write_csv(final_symptom_onsets_share, path="augmented_data/augmented_symptom_times.csv")
  write_csv(final_infections_share, path="augmented_data/augmented_infection_times.csv")
  
  ## Need to free some memory
  rm(final_dat_to_share)
  rm(final_symptom_onsets_share)
  rm(final_infections_share)
  
  final_dat_to_share_province <- final_all_province %>% filter(repeat_no %in% use_repeats)
  colnames(final_dat_to_share_province) <- c("repeat_no","variable","date","province","from_confirmed","inflated_event","total_events")
  
  write_csv(final_dat_to_share_province, path="augmented_data/augmented_totals_by_province.csv")
}