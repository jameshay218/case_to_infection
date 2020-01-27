setwd("~/GitHub/case_to_infection/")

library(ggplot2)
library(tidyverse)
library(lubridate)
library(cowplot)
library(tibble)
library(patchwork)
library(ggpmisc)
library(ggpubr)
library(plyr)
source("code/analysis_functions.R")

hubei_data_path <- "~/GitHub/nCoV2019/ncov_hubei.csv"
other_data_path <- "~/GitHub/nCoV2019/ncov_outside_hubei.csv"

key_colnames <- c("ID","age","country", "sex","city","province",
                  "latitude","longitude")
var_colnames <- c("date_confirmation","date_onset_symptoms","date_admission_hospital")
use_colnames <- c(key_colnames, var_colnames)

mean_incubation <- 10
var_incubation <- 12.33
repeats <- 100

## First step is to clean and take a look at the data
## This combines the data for Hubei and other locations in China
## Note this is ONLY China
source("code/data_clean_all.R")
#p_other_data
#p_other_hosp
#p_other_confirm

####################################
## CONFIRMATION DELAY DISTRIBUTION
####################################
## Assume that there is at least a 1 day delay to reporting, so < 1 day is set to 1
other_dat2[other_dat2$confirmation_delay < 1, "confirmation_delay"] <- 1

## Fit a geometric distribution to the confirmation delay distribution
fit1 <- optim(c(0.1), fit_geometric, dat=other_dat2$confirmation_delay-1,method="Brent",lower=0,upper=1)
fit_line1 <- dgeom(seq(0,25,by=1),prob=fit1$par)
fit_line_dat1 <- data.frame(x=seq(1,26,by=1),y=fit_line1)

p_other_confirm_fit<- ggplot(other_dat2) + 
  geom_histogram(aes(x=confirmation_delay,y=..density..),binwidth=1,col="black") +
  geom_line(data=fit_line_dat1, aes(x=x,y=y), col="red",size=1) +
  scale_x_continuous(breaks=seq(0,30,by=5),labels=seq(0,30,by=5)) +
  scale_y_continuous(expand=c(0,0),limits=c(0,0.2)) +
  geom_vline(xintercept=1,linetype="dashed") +
  ylab("Probability density") + xlab("Days since symptom onset") +
  ggtitle("Distribution of delays between symptom\n onset and confirmation") +
  theme_pubr()
#p_other_confirm_fit
####################################
## HOSPITALISATION DELAY DISTRIBUTION
####################################
## Fit a geometric distribution to hospitalisation delay distribution
fit2 <- optim(c(0.1), fit_geometric, dat=other_dat1$hospitalisation_delay,method="Brent",lower=0,upper=1)
times <- seq(0,25,by=1)
fit_line2 <- dgeom(times,prob=fit2$par)
fit_line_dat2 <- data.frame(x=times,y=fit_line2)

p_other_hosp_fit<- ggplot(other_dat1) + 
  geom_histogram(aes(x=hospitalisation_delay,y=..density..),binwidth=1) +
  geom_line(data=fit_line_dat2, aes(x=x,y=y), col="red") +
  scale_x_continuous(breaks=seq(0,25,by=1)) +
  ggtitle("Distribution of delays between symptom\n onset and hospitalisation (not great fit)") +
  theme_pubr()
## Fit isn't great for first day

####################################
## SYMPTOM ONSET DISTRIBUTION
####################################
## Going to use a gamma distibution with mode of 10 days
times <- seq(0,25,by=0.1)
incubation_period <- dgamma_mean(times,mean_incubation, var_incubation, FALSE)
inc_data <- data.frame(day=times,probability=incubation_period)
p_incubation <- ggplot(inc_data) + 
  geom_ribbon(aes(x=times, ymax=incubation_period,ymin=0),fill="grey70",col="black",size=1) + 
  ylab("Probability density") +
  xlab("Days since onset of infection") +
  ggtitle("Incubation period distribution\n (time from infection to symptoms)") +
  scale_y_continuous(limits=c(0,0.15),expand=c(0,0),breaks=seq(0,0.15,by=0.05)) +
  scale_x_continuous(expand=c(0,0)) +
  theme_pubr()

#assumption_plot
## Now for each reported case without symptom onset time, going to generate a random
## symptom onset time from the geometric distribution
other_dat_china <- other_dat[other_dat$country == "China",]
other_dat_china$augmented <- "No"
which_to_sim <- which(is.na(other_dat_china$date_onset_symptoms) & !is.na(other_dat_china$date_confirmation))
other_dat_china[which_to_sim,"augmented"] <- "Augmented symptom onset and infection"
n_to_sim <- nrow(other_dat_china[which_to_sim,])

sim_confirmation_delays <- rgeom(n_to_sim, fit1$par) + 1
other_dat_china[which_to_sim,"date_onset_symptoms"] <- other_dat_china[which_to_sim,"date_confirmation"] - sim_confirmation_delays

which_to_sim_infection <- which(!is.na(other_dat_china$date_onset_symptoms))
sim_incubation_times <- rgamma_mean(nrow(other_dat_china[which_to_sim_infection,]), mean_incubation, var_incubation)
other_dat_china[setdiff(which_to_sim_infection, which_to_sim),"augmented"] <- "Augmented infection only"
other_dat_china$date_infection <- other_dat_china$date_onset_symptoms
other_dat_china[which_to_sim_infection,"date_infection"] <- other_dat_china[which_to_sim_infection,"date_onset_symptoms"] - sim_incubation_times
other_dat_china$augmented <- as.factor(other_dat_china$augmented)
other_dat_tmp2 <- reshape2::melt(other_dat_china[,c("age","sex","country","hubei",
                                              "date_onset_symptoms","date_admission_hospital",
                                              "date_infection","date_confirmation","augmented")], 
                                 id.vars=c("age","sex","country","hubei","augmented"))
other_dat_tmp2 <- other_dat_tmp2[!is.na(other_dat_tmp2$value),]
other_dat_tmp2[other_dat_tmp2$variable == "date_confirmation","augmented"] <- "No"
other_dat_tmp2[other_dat_tmp2$variable == "hospital_admission_date","augmented"] <- "No"

variable_key <- c("date_confirmation"="Confirmation date (known)",
                  "date_onset_symptoms"="Onset of symptoms",
                  "date_admission_hospital"="Hospital admission date",
                  "date_infection"="Augmented infection date")
other_dat_tmp2$variable <- as.character(other_dat_tmp2$variable)
other_dat_tmp2$variable <- variable_key[other_dat_tmp2$variable]

p_data_augmented_example <- ggplot(other_dat_tmp2) + 
  geom_histogram(aes(x=value, fill=augmented),binwidth=1,col="black") +
  scale_fill_manual(values=c("No"="grey40","Augmented infection only"="orange","Augmented symptom onset and infection"="blue")) +
  xlab("Date") + ylab("Count") +
  ggtitle("Single simulation of augmented infection time data") +
  facet_wrap(~variable, ncol=2) +
  theme_bw() +
  scale_x_date(limits=c(convert_date("01.12.2019"),convert_date(latest_date)),
               breaks="7 day") +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.position = "bottom")
p_data_augmented_example

## FULL AUGMENTATION
other_dat_china1 <- other_dat[other_dat$country == "China",]

## Now let's repeat this process many times to get a distribution
sim_data_infections <- matrix(NA, nrow=repeats, ncol=nrow(other_dat_china))
sim_data_symptoms <- matrix(NA, nrow=repeats, ncol=nrow(other_dat_china))

for(i in seq_len(repeats)){
  tmp <- augment_infection_times(other_dat_china1, mean_incubation, var_incubation, fit1$par)
  sim_data_infections[i,] <- tmp$augmented_infection_times
  sim_data_symptoms[i,] <- tmp$augmented_symptom_onsets
}

sim_data_infections_melted <- reshape2::melt(sim_data_infections)
sim_data_symptoms_melted <- reshape2::melt(sim_data_symptoms)

sim_data_infections_melted$var <- "date_infection"
sim_data_symptoms_melted$var <- "date_onset_symptoms"

colnames(sim_data_infections_melted) <- colnames(sim_data_symptoms_melted) <- c("repeat_no","individual","date","var")
sim_data_all <- rbind(sim_data_infections_melted, sim_data_symptoms_melted)
sim_data_all$date <- as.Date(floor(sim_data_all$date), origin="1970-01-01")

sim_data_sum <- ddply(sim_data_all, .(repeat_no, var, date), nrow)

variable_key2 <- c("date_confirmation"="Confirmation date (known)",
                   "date_onset_symptoms"="Onset of symptoms for cases observed to date",
                   "date_admission_hospital"="Hospital admission date",
                   "date_infection"="Augmented infection date for cases observed to date")

sim_data_zeros <- expand.grid(repeat_no=unique(sim_data_sum$repeat_no),
                              var=unique(sim_data_sum$var),
                              date=unique(sim_data_sum$date))
sim_data_sum <- merge(sim_data_sum, sim_data_zeros,all=TRUE)
sim_data_sum[is.na(sim_data_sum$V1),"V1"] <- 0

################################################
## OVERALL PLOT
## Distribution of times for each date
sim_data_quantiles <- ddply(sim_data_sum, .(date, var), function(x) quantile(x$V1, c(0.025,0.5,0.975),na.rm=TRUE))

## Get confirmation time data
confirm_data <- ddply(other_dat_china1[!is.na(other_dat_china1$date_confirmation),], ~date_confirmation, function(x) nrow(x))
confirm_data$Variable <- "Confirmed cases"

sim_data_quantiles$var <- variable_key2[sim_data_quantiles$var]

colnames(sim_data_quantiles) <- c("date","Variable","lower","median","upper")
augmented_data_plot <- plot_augmented_data(sim_data_quantiles, confirm_data)

#assumption_plot
#augmented_data_plot
#pdf("fig/augmented_data_plot.pdf",height=8,width=10)
#augmented_data_plot
#dev.off()


#pdf("fig/assumption_plot.pdf",height=5,width=12)
#assumption_plot
#dev.off()
## Distribution of times for each individual
sim_data_quantiles_indiv <- ddply(sim_data_all, .(individual, var), function(x) quantile(as.numeric(x$date), c(0.025,0.5,0.975),na.rm=TRUE))
sim_data_quantiles_indiv$`2.5%` <- as.Date(sim_data_quantiles_indiv$`50%`, origin="1970-01-01")
sim_data_quantiles_indiv$`50%` <- as.Date(sim_data_quantiles_indiv$`2.5%`, origin="1970-01-01")
sim_data_quantiles_indiv$`97.5%` <- as.Date(sim_data_quantiles_indiv$`97.5%`, origin="1970-01-01")



#######################
## SPATIAL PLOTS
#######################
individual_key <- other_dat_china1[,c("ID","age","country","sex","city","province","latitude","longitude")]
colnames(individual_key)[1] <- "individual"

#sim_data_all_wide <- sim_data_all %>% 
#  pivot_wider(id_cols=c("repeat_no","individual"),names_from="var",
#              values_from="date")
merged_data <- merge(individual_key, sim_data_all)
merged_data <- merged_data[!is.na(merged_data$date),]
#############################
## Aggregate by province
## Get confirmation time data
confirm_data_province <- ddply(other_dat_china1[!is.na(other_dat_china1$date_confirmation),], .(province,date_confirmation), function(x) nrow(x))
confirm_data_province$Variable <- "Confirmed cases"
province_data <- ddply(merged_data, .(repeat_no, var, date, province), nrow)
sim_data_zeros_prov <- expand.grid(repeat_no=unique(province_data$repeat_no),
                              var=unique(province_data$var),
                              province=unique(province_data$province),
                              date=unique(province_data$date))
province_data <- merge(province_data, sim_data_zeros_prov,all=TRUE)
province_data[is.na(province_data$V1),"V1"] <- 0

sim_data_quantiles_province <- ddply(province_data, .(date, var, province), 
                                     function(x) quantile(x$V1, c(0.025,0.5,0.975),na.rm=TRUE))

sim_data_quantiles_province$var <- variable_key2[sim_data_quantiles_province$var]
colnames(sim_data_quantiles_province) <- c("date","Variable","province","lower","median","upper")

total_confirmed_prov <- ddply(confirm_data_province, ~province, function(x) sum(x$V1))
total_confirmed_prov <- total_confirmed_prov[order(-total_confirmed_prov$V1),]
factor_order <- as.character(total_confirmed_prov$province)

confirm_data_province$province <- factor(as.character(confirm_data_province$province), 
                                            levels=factor_order)
sim_data_quantiles_province$province <- factor(as.character(sim_data_quantiles_province$province), 
                                         levels=factor_order)



by_province <- plot_augmented_data_province(sim_data_quantiles_province, confirm_data_province)
top_6 <- factor_order[1:6]
by_province_top6 <- plot_augmented_data_province(sim_data_quantiles_province[sim_data_quantiles_province$province %in% top_6,], 
                                                 confirm_data_province[confirm_data_province$province %in% top_6,])
by_province_top6 <- by_province_top6 + facet_wrap(~province, ncol=3, scales="free_y") + theme(legend.text=element_text(size=10))
by_province_top6
#########################
## FINAL HOUSEKEEPING
## Tidy up data to share
sim_data_infections1 <- as.data.frame(t(sim_data_infections[1:100,]))
for(i in seq_len(ncol(sim_data_infections1))){
  sim_data_infections1[,i] <- as.Date(floor(sim_data_infections1[,i]), origin="1970-01-01")
}
sim_data_infections1 <- cbind(other_dat_china1, sim_data_infections1)


sim_data_symptoms1 <- as.data.frame(t(sim_data_symptoms[1:100,]))
for(i in seq_len(ncol(sim_data_symptoms1))){
  sim_data_symptoms1[,i] <- as.Date(floor(sim_data_symptoms1[,i]), origin="1970-01-01")
}
sim_data_symptoms1 <- cbind(other_dat_china1, sim_data_symptoms1)

write.csv(sim_data_infections1, "augmented_infection_times.csv")
write.csv(sim_data_symptoms1, "augmented_symptom_times.csv")


## Create results panel plot programmatically
element_text_size <- 10
text_size_theme <- theme(title=element_text(size=element_text_size), 
                         axis.text=element_text(size=element_text_size), 
                         axis.title = element_text(size=element_text_size))
p_other_confirm_fit1 <- p_other_confirm_fit + text_size_theme
p_incubation1 <- p_incubation + text_size_theme
assumption_plot <- plot_grid(p_other_confirm_fit1, p_incubation1,ncol=2,align="hv")
augmented_data_plot1 <- augmented_data_plot + theme(legend.position=c(0.25,0.25))
layout <- c(
  area(t=0,b=12,l=0,r=18),
  area(t=2,b=7,l=2,r=14)
)

results_panel <- augmented_data_plot1 + assumption_plot + plot_layout(design=layout)

#pdf("fig/results_plot.pdf", height=8, width=12)
#results_panel
#dev.off()

#png("fig/results_plot.png", height=8, width=12, units="in",res=300)
#results_panel
#dev.off()

#pdf("fig/by_province.pdf", height=12, width=12)
#by_province
#dev.off()

#png("fig/by_province.png", height=12, width=12, units="in",res=300)
#by_province
#dev.off()

#png("fig/by_province_top6.png", height=8, width=10, units="in",res=300)
#by_province_top6
#dev.off()


