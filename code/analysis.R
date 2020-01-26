setwd("~/GitHub/case_to_infection/")

library(ggplot2)
library(tidyverse)
library(lubridate)
library(cowplot)
library(ggpubr)
library(plyr)
source("code/analysis_functions.R")

## Assumed parameters
mean_incubation <- 10
var_incubation <- 12.33
repeats <- 1000

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
  scale_x_continuous(breaks=seq(0,5,by=1)) +
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

assumption_plot <- plot_grid(p_other_confirm_fit, p_incubation,ncol=2,align="hv")
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
  scale_x_date(limits=c(convert_date("01.12.2019"),convert_date("31.01.2020")),
               breaks="7 day") +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.position = "bottom")
#p_data_augmented_example

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


## Distribution of times for each date
sim_data_quantiles <- ddply(sim_data_sum, .(date, var), function(x) quantile(x$V1, c(0.025,0.5,0.975),na.rm=TRUE))

## Get confirmation time data
confirm_data <- ddply(other_dat_china1[!is.na(other_dat_china1$date_confirmation),], ~date_confirmation, function(x) nrow(x))
confirm_data$Variable <- "Confirmed cases"

variable_key2 <- c("date_confirmation"="Confirmation date (known)",
                  "date_onset_symptoms"="Onset of symptoms for cases observed to date",
                  "date_admission_hospital"="Hospital admission date",
                  "date_infection"="Augmented infection date for cases observed to date")

sim_data_quantiles$var <- variable_key2[sim_data_quantiles$var]


colnames(sim_data_quantiles) <- c("date","Variable","lower","median","upper")
augmented_data_plot <- ggplot(sim_data_quantiles) +
  geom_bar(data=confirm_data,aes(x=date_confirmation,y=V1,fill=Variable),stat="identity") +
  geom_ribbon(aes(x=date,ymax=upper,ymin=lower,fill=Variable,col=Variable),alpha=0.25) +
  geom_line(aes(x=date, y=median,col=Variable),size=1) +
  scale_y_continuous(limits=c(0,300),expand=c(0,0),breaks=seq(0,300,by=25)) +
  scale_x_date(limits=c(convert_date("01.12.2019"),convert_date("26.01.2020")),
               breaks="5 day") + 
  scale_fill_manual(values=c("orange","grey40","blue")) + scale_color_manual(values=c("orange","blue")) +
  ggtitle("Augmented and observed timings of infection and symptom onset in China") +
  ylab("Count") + xlab("Date of event") +
  theme_pubr() +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        panel.grid.major = element_line(colour="grey70"),
        legend.position = "none") 
#assumption_plot
#augmented_data_plot
pdf("augmented_data_plot.pdf",height=8,width=10)
augmented_data_plot
dev.off()


pdf("assumption_plot.pdf",height=5,width=12)
assumption_plot
dev.off()
## Distribution of times for each individual
sim_data_quantiles_indiv <- ddply(sim_data_all, .(individual, var), function(x) quantile(as.numeric(x$date), c(0.025,0.5,0.975),na.rm=TRUE))
sim_data_quantiles_indiv$`2.5%` <- as.Date(sim_data_quantiles_indiv$`50%`, origin="1970-01-01")
sim_data_quantiles_indiv$`50%` <- as.Date(sim_data_quantiles_indiv$`2.5%`, origin="1970-01-01")
sim_data_quantiles_indiv$`97.5%` <- as.Date(sim_data_quantiles_indiv$`97.5%`, origin="1970-01-01")



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

