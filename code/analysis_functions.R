## Vector of confounding date entries
confounding_dates <- c("pre 18.01.2020",
                       "early january",
                       "",
                       "none")

## Checks the given vector of dates as strings 
## and converts to NA if not useable
clean_dates <- function(dates){
  dates_new <- dates
  dates_new[dates_new %in% confounding_dates] <- NA
  dates_new
  
}

convert_date <- function(date, format="%d.%m.%Y", origin="01.01.1970"){
  as.Date(date, origin=origin, tryFormat=c(format))
}


fit_geometric <- function(prob, dat){
  -sum(dgeom(x=dat, prob, log=TRUE))
}

dgamma_mean <- function(x, mean, var, use_log=FALSE){
  scale <- var/mean
  shape <- mean/scale
  probs <- dgamma(x,shape=shape,scale=scale,log=use_log)
  probs
}

rgamma_mean <- function(n, mean, var){
  scale <- var/mean
  shape <- mean/scale
  rgamma(n,shape=shape,scale=scale)
}

fit_gamma <- function(pars, dat){
  mean <- pars[1]
  var <- pars[2]
  sum(dgamma_mean(dat, mean, var, TRUE))
}


augment_infection_times <- function(dat, mean_incubation, var_incubation, geom_prob){
  which_to_sim <- which(is.na(dat$date_onset_symptoms) & !is.na(dat$date_confirmation))
  n_to_sim <- nrow(dat[which_to_sim,])
  sim_confirmation_delays <- rgeom(n_to_sim, geom_prob) + 1
  dat[which_to_sim,"date_onset_symptoms"] <- dat[which_to_sim,"date_confirmation"] - sim_confirmation_delays
  
  
  which_to_sim_infection <- which(!is.na(other_dat_china$date_onset_symptoms))
  sim_incubation_times <- rgamma_mean(nrow(other_dat_china[which_to_sim_infection,]), mean_incubation, var_incubation)
  dat$date_infection <- dat$date_onset_symptoms
  dat[which_to_sim_infection,"date_infection"] <- dat[which_to_sim_infection,"date_onset_symptoms"] - sim_incubation_times
  
  return(list(augmented_symptom_onsets=dat$date_onset_symptoms,
              augmented_infection_times=dat$date_infection   
         ))
  
  
  
}
