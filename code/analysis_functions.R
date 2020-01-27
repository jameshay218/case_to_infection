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
  dat[which_to_sim,"date_onset_symptoms"] <- dat[which_to_sim,"date_confirmation"] - floor(sim_confirmation_delays)
  
  
  which_to_sim_infection <- which(!is.na(other_dat_china$date_onset_symptoms))
  sim_incubation_times <- rgamma_mean(nrow(other_dat_china[which_to_sim_infection,]), mean_incubation, var_incubation)
  dat$date_infection <- dat$date_onset_symptoms
  dat[which_to_sim_infection,"date_infection"] <- dat[which_to_sim_infection,"date_onset_symptoms"] - floor(sim_incubation_times)
  
  return(list(augmented_symptom_onsets=dat$date_onset_symptoms,
              augmented_infection_times=dat$date_infection   
         ))
}


plot_augmented_data <- function(data_quantiles, confirmed_data, max_date="27.01.2020"){
  p <- ggplot(data_quantiles) +
    geom_bar(data=confirmed_data,aes(x=date_confirmation,y=V1,fill=Variable),stat="identity") +
    geom_ribbon(aes(x=date,ymax=upper,ymin=lower,fill=Variable,col=Variable),alpha=0.25) +
    geom_line(aes(x=date, y=median,col=Variable),size=1) +
    scale_y_continuous(limits=c(0,300),expand=c(0,0),breaks=seq(0,300,by=25)) +
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
