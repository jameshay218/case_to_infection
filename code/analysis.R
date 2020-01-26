setwd("~/GitHub/case_to_infection/")

library(ggplot2)
library(tidyverse)
library(lubridate)
library(cowplot)
source("code/analysis_functions.R")

## First step is to take a look at the data
source("code/data_clean_outside_hubei.R")
p_other_data
p_other_hosp
p_other_confirm

source("code/data_clean_hubei.R")
p_hubei_data
p_hubei_dist

####################################
## CONFIRMATION DELAY DISTRIBUTION
####################################

## Will use the outside hubei data, as much more complete
## Big assumption that confirmation delays are similar!

## Assume that there is at least a 1 day delay to reporting, so < 1 day is set to 1
other_dat2[other_dat2$confirmation_delay < 1, "confirmation_delay"] <- 1

## Fit a geometric distribution to this
fit1 <- optim(c(0.1), fit_geometric, dat=other_dat2$confirmation_delay-1,method="Brent",lower=0,upper=1)
fit_line1 <- dgeom(seq(0,30,by=1),prob=fit1$par)
fit_line_dat1 <- data.frame(x=seq(1,31,by=1),y=fit_line1)

p_other_confirm_fit<- ggplot(other_dat2) + 
  geom_histogram(aes(x=confirmation_delay,y=..density..),binwidth=1) +
  geom_line(data=fit_line_dat1, aes(x=x,y=y), col="red") +
  scale_x_continuous(breaks=seq(0,max(other_dat2$confirmation_delay + 10),by=1)) +
  ggtitle("Distribution of delays between symptom onset and confirmation") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        panel.grid.minor=element_blank())
p_other_confirm_fit


####################################
## HOSPITALISATION DELAY DISTRIBUTION
####################################
## Will use the outside hubei data, as much more complete
## Big assumption that hospitalisation delays are similar!

## Fit a geometric distribution to this
fit2 <- optim(c(0.1), fit_geometric, dat=other_dat1$hospitalisation_delay,method="Brent",lower=0,upper=1)
times <- seq(0,15,by=1)
fit_line2 <- dgeom(times,prob=fit2$par)
fit_line_dat2 <- data.frame(x=times,y=fit_line2)

p_other_hosp_fit<- ggplot(other_dat1) + 
  geom_histogram(aes(x=hospitalisation_delay,y=..density..),binwidth=1) +
  geom_line(data=fit_line_dat2, aes(x=x,y=y), col="red") +
  scale_x_continuous(breaks=seq(0,max(other_dat1$hospitalisation_delay + 10),by=1)) +
  ggtitle("Distribution of delays between symptom onset and hospitalisation (not great fit)") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        panel.grid.minor=element_blank())
## Fit isn't great for first day
p_other_hosp_fit



####################################
## SYMPTOM ONSET DISTRIBUTION
####################################
## Going to use a gamma distibution with mode of 10 days
mean_incubation <- 10
var_incubation <- 2
times <- seq(0,25,by=0.1)
incubation_period <- dgamma_mean(times,mean_incubation, var_incubation, FALSE)
plot(incubation_period~times, xlab="Day of symptom onset from infection",ylab="Density", type="l")


## Now for each reported case without symptom onset time, going to generate a random
## symptom onset time from the geometric distribution
hubei_dat[is.na(hubei_dat$date_onset_symptoms) & !is.na(hubei_dat$date_admission_hospital),]





