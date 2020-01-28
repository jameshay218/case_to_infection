#' Vector of confounding date entries
confounding_dates <- c("none", "10.01.2020 - 22.01.2020", "pre 18.01.2020", "early january","",
  "11.26.2020", "18.01.2020 - 23.01.2020", "not sure")

#' Checks the given vector of dates as strings and converts to NA if not useable
clean_dates <- function(dates){
  dates_new <- dates
  dates_new[dates_new %in% confounding_dates] <- NA
  dates_new
  
}
#' Converte date string to date format used in analysis
convert_date <- function(date1, format="%d.%m.%Y", origin="1970-01-01"){
  if(!(all(is.numeric(date1[!is.na(date1)])))) {
      ## Clean provided vector of dates
    date1 <- sapply(date1, function(x){
      clean <- x
      ## Keep as NA if needed
      if(!is.na(x)){
        ## Check for format dd.mm.yyyy or dd/mm/yyyy etc
        if(nchar(x) > 10){
          message(cat("Provided date \"", x, "\" must be less than 11 characters long\n"))
          clean <- NA
        } else {
          ## See if we can convert date. If it's something like "none", then will return NA
          clean <- tryCatch(as.Date(x, origin=origin, tryFormat=c(format)),
                                    error = function(e) {
                                      message(cat("Could not convert date \"", x, "\"\n"))
                                      NA})
        }
      }
      clean
    }, simplify=FALSE)
    ## Seems redundant, but converting back to a vector needs going via numeric vector
    date1 <- as.numeric(unlist(date1))
  }
  as.Date(date1, origin=origin)
}

#' Function for optim to fit a geometric distribution
#' 
#' @param prob probability of success
#' @param dat is a vector of event observation times
#' @return negative sum log likelihood from dgeom
fit_geometric <- function(prob, dat){
  -sum(dgeom(x=dat, prob, log=TRUE))
}

#' Gamma prob density from mean and variance
dgamma_mean <- function(x, mean, var, use_log=FALSE){
  scale <- var/mean
  shape <- mean/scale
  probs <- dgamma(x,shape=shape,scale=scale,log=use_log)
  probs
}

#' Gamma random generator from mean and variance
rgamma_mean <- function(n, mean, var){
  scale <- var/mean
  shape <- mean/scale
  rgamma(n,shape=shape,scale=scale)
}
#' Function for optim to fit a gamma distribution
#' 
#' @param pars vector, 1: gamma mean; 2: gamma variance
#' @param dat is a vector of event observation times
#' @return negative sum log likelihood from dgamma_mean
fit_gamma <- function(pars, dat){
  mean <- pars[1]
  var <- pars[2]
  sum(dgamma_mean(dat, mean, var, TRUE))
}

#' Augment infection times
augment_infection_times <- function(dat, inc_period_alpha, inc_period_sigma, p_confirm_delay){
  ## Sim symptom onset times for those without
  which_to_sim <- which(is.na(dat$date_onset_symptoms) & !is.na(dat$date_confirmation))
  n_to_sim <- nrow(dat[which_to_sim,])
  
  ## From the provided geometric distribution + 1
  sim_confirmation_delays <- rgeom(n_to_sim, p_confirm_delay) + 1
  dat[which_to_sim,"date_onset_symptoms"] <- dat[which_to_sim,"date_confirmation"] - floor(sim_confirmation_delays)
  
  ## Sim infection times for everyone with symptom onset time
  which_to_sim_infection <- which(!is.na(dat$date_onset_symptoms))
  ## Old 
  ## sim_incubation_times <- rgamma_mean(nrow(other_dat_china[which_to_sim_infection,]), mean_incubation, var_incubation)
  
  ## New
  sim_incubation_times <- rweibull(nrow(dat[which_to_sim_infection,]), alpha, sigma)
  dat$date_infection <- dat$date_onset_symptoms
  dat[which_to_sim_infection,"date_infection"] <- dat[which_to_sim_infection,"date_onset_symptoms"] - floor(sim_incubation_times)
  
  return(list(augmented_symptom_onsets=dat$date_onset_symptoms,
              augmented_infection_times=dat$date_infection   
         ))
}


plot_augmented_data <- function(data_quantiles, confirmed_data, max_date="27.01.2020",
                                ymax=500,ybreaks=25){
  p <- ggplot(data_quantiles) +
    geom_bar(data=confirmed_data,aes(x=date_confirmation,y=V1,fill=Variable),stat="identity") +
    geom_ribbon(aes(x=date,ymax=upper,ymin=lower,fill=Variable,col=Variable),alpha=0.25) +
    geom_line(aes(x=date, y=median,col=Variable),size=1) +
    scale_y_continuous(limits=c(0,ymax),expand=c(0,0),breaks=seq(0,ymax,by=ybreaks)) +
    scale_x_date(limits=c(convert_date("01.12.2019"),convert_date(max_date)),
                 breaks="5 day") + 
    scale_fill_manual(values=c("orange","grey40","blue")) + 
    scale_color_manual(values=c("orange","blue"),guide="none") +
    ggtitle("Augmented and observed timings of infection and symptom onset in China") +
    ylab("Count") + xlab("Date of event") +
    theme_pubr() +
    theme(axis.text.x=element_text(angle=45,hjust=1),
          panel.grid.major = element_line(colour="grey70"),
          legend.position = c(0.25,0.75)) 
  p
}

plot_augmented_data_province <- function(data_quantiles_province, confirmed_data_province, max_date="27.01.2020"){
  p <- ggplot(data_quantiles_province) +
  geom_bar(data=confirmed_data_province,aes(x=date_confirmation,y=V1,fill=Variable),stat="identity") +
  geom_ribbon(aes(x=date,ymax=upper,ymin=lower,fill=Variable,col=Variable),alpha=0.25) +
  geom_line(aes(x=date, y=median,col=Variable),size=1) +
  scale_x_date(limits=c(convert_date("01.12.2019"),convert_date(max_date)),
               breaks="7 day") + 
  scale_fill_manual(values=c("orange","grey40","blue")) + 
  scale_color_manual(values=c("orange","blue"),guide="none") +
  ggtitle("Augmented and observed timings of infection and symptom onset in China by province\n ordered by total confirmed cases") +
    geom_hline(yintercept=0,linetype="dashed",col="grey80",size=0.5) +
  ylab("Count") + xlab("Date of event") +
  facet_wrap(~province, scales="free_y", ncol=5) +
  theme_pubr() +
  theme(axis.text.x=element_text(angle=45,hjust=1,size=12),
        axis.text.y=element_text(size=12),
        title=element_text(size=14),
        strip.text = element_text(size=12),
        legend.text=element_text(size=12),
        legend.position = "bottom") 
  p
}
