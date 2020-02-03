#' Vector of confounding date entries
confounding_dates <- c("none", "10.01.2020 - 22.01.2020", "pre 18.01.2020", "early january","",
  "11.26.2020", "18.01.2020 - 23.01.2020", "not sure")

#' Colors for plots
blue_color <- "#0072B2"
orange_color <- "#D55E00"

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
  -sum(dgamma_mean(dat, mean, var, TRUE))
}

#' Function for optim to fit poisson distribution
fit_poisson <- function(pars, dat){
  mean <- pars[1]
  -sum(dpois(dat, mean, TRUE))
}

dgamma_discrete_mean <- function(dat, mean, var, use_log=TRUE){
  scale <- var/mean
  shape <- mean/scale
  ddgamma(dat, shape=shape,scale=scale,log=use_log)
}

#' Discretised gamma
dgamma_mean_discrete <- function(dat, mean, var, use_log=TRUE){
  dat_low <- which(dat <= 1)
  scale <- var/mean
  shape <- mean/scale
  
  y <- numeric(length(dat))
  if(use_log){
    y <- log(pgamma(dat+1, shape=shape, scale=scale, lower = TRUE) - pgamma(dat, shape=shape, scale=scale, lower = TRUE))
    y[dat_low] <- pgamma(1,  shape=shape, scale=scale,  lower = TRUE,log=TRUE)
  } else {
    y <- pgamma(dat+1,  shape=shape, scale=scale,  lower = TRUE) - pgamma(dat, shape=shape, scale=scale, lower = TRUE)
    y[dat_low] <- pgamma(1,  shape=shape, scale=scale,  lower = TRUE)
  }
  y
  
}
fit_gamma_discrete <- function(pars,dat){
  mean <- pars[1]
  var <- pars[2]
  scale <- var/mean
  shape <- mean/scale
  -sum(ddgamma(dat, shape=shape, scale=scale, log=TRUE))
  
}

#' Augment infection times
augment_infection_times <- function(dat, inc_period_alpha, inc_period_sigma, gamma_mean, gamma_var){
  ## Sim symptom onset times for those without
  which_to_sim <- which(is.na(dat$date_onset_symptoms) & !is.na(dat$date_confirmation))
  n_to_sim <- nrow(dat[which_to_sim,])
  
  ## From the provided geometric distribution + 1
  #sim_confirmation_delays <- rgeom(n_to_sim, p_confirm_delay) + 1
  sim_confirmation_delays <- rgamma_mean(n_to_sim, gamma_mean, gamma_var) 
  dat[which_to_sim,"date_onset_symptoms"] <- dat[which_to_sim,"date_confirmation"] - ceiling(sim_confirmation_delays)
  
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

#' Plot time from first infection in Hubei by province
#' 
plot_time_from_start <- function(sim_data_infections_melted, individual_key,xmax=100) {
  melted_infs <- sim_data_infections_melted %>% as_tibble() %>% select(-var)
  
  melted_infs <- right_join(melted_infs, individual_key,by="individual")
  start_inf_hubei <- melted_infs %>% group_by(repeat_no) %>% filter(province=="Hubei") %>% summarise(hubei_time=min(date,na.rm=TRUE))
  start_times <- melted_infs %>% group_by(repeat_no, province) %>% summarise(start_time=min(date,na.rm=TRUE))
  
  all_data <- left_join(start_inf_hubei, start_times,by="repeat_no")
  all_data <- all_data %>% mutate(delay=start_time-hubei_time)
  
  ## Order by delay
  orders <- all_data %>% group_by(province) %>% summarise(mean_start=mean(delay)) %>% arrange(-mean_start) %>% pull(province)
  all_data$province <- factor(all_data$province, levels=orders)
  ggplot(all_data) + 
    geom_violin(aes(y=delay, x=province),
                draw_quantiles=c(0.025,0.5,0.975),scale="width",trim = TRUE,
                fill="grey80") +
    coord_flip() +
    xlab("Province") + ylab("Delay between first infection in Hubei \nand first infection in province (days)") +
    scale_y_continuous(limits=c(0,xmax)) +
    theme_bw() 
  
}


plot_augmented_data <- function(data_quantiles, confirmed_data, max_date="27.01.2020",
                                min_date="01.12.2019", cols = c("grey40","blue","skyblue","orange","red"),
                                cols2 = c("blue","skyblue","orange","red"),
                                ymax=500,ybreaks=25, thresholds=NULL,threshold_y=700,
                                title = "Augmented and observed timings of infection and symptom onset in China"){
  threshold_dat <- data.frame(xmin=c(convert_date("01.12.2019"), thresholds),
                              xmax=c(thresholds, convert_date(max_date)),
                              fills=c(">95%",">80%",">50%",">20%","<20%"))

  p <- ggplot(data_quantiles)
  if(!is.null(thresholds)) {
    threshold_dat <- data.frame(text=c(">95%",">80%",">50%",">20%"),
                                x_val=thresholds)
    p <- p +    
      geom_vline(xintercept=thresholds,linetype="dashed") +
      geom_label(data=threshold_dat,aes(x=x_val, y=ymax*0.4, label=text))
  }
  p <- p +
    #geom_rect(data=threshold_dat,aes(xmin=xmin,xmax=xmax,ymin=0,ymax=ymax,alpha=fills),fill="red") +
    geom_bar(data=confirmed_data,aes(x=date_confirmation,y=n,fill=Variable),stat="identity") +
    geom_ribbon(aes(x=date,ymax=upper,ymin=lower,fill=Variable,col=Variable),alpha=0.25) +
    geom_line(aes(x=date, y=mean,col=Variable),size=1) +
    scale_y_continuous(expand=c(0,0),breaks=seq(0,ymax,by=ybreaks)) +
    coord_cartesian(ylim=c(0,ymax),xlim=c(convert_date(min_date), convert_date(max_date)+1)) +
    scale_x_date(limits=c(convert_date("01.12.2019"),convert_date(max_date)+1),
                 breaks="5 day") + 
    scale_fill_manual(values= cols) + 
    scale_color_manual(values=cols2,guide="none") +
    ggtitle(title) +
    ylab("Count") + xlab("Date of event") +
    theme_pubr() +
    theme(axis.text.x=element_text(angle=45,hjust=1),
          panel.grid.major = element_line(colour="grey70"),
          legend.position = c(0.25,0.75)) 
  p
}


plot_augmented_data_province <- function(data_quantiles_province, confirmed_data_province, max_date="27.01.2020"){
  p <- ggplot(data_quantiles_province) +
  geom_bar(data=confirmed_data_province,aes(x=date_confirmation,y=n,fill=Variable),stat="identity") +
  geom_ribbon(aes(x=date,ymax=upper,ymin=lower,fill=Variable,col=Variable),alpha=0.25) +
  geom_line(aes(x=date, y=median,col=Variable),size=1) +
  scale_x_date(limits=c(convert_date("01.12.2019"),convert_date(max_date)),
               breaks="7 day") + 
  scale_fill_manual(values=c("blue","grey40","orange")) + 
  scale_color_manual(values=c("blue","orange"),guide="none") +
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

expand_arcgis <- function(arcgis_dat) {
  arcgis_dat %>% group_by(date_confirmation, province, country) %>% filter(n > 0) %>% drop_na() %>% uncount(n)
}

merge_data <- function(linelist_dat, arcgis_dat, switch_date){
  arcgis_dat <- arcgis_dat %>% select(province, raw_day, country_region, country, diff)
  colnames(arcgis_dat) <- c("province","date_confirmation","country_region","country","n")
  arcgis_dat <- arcgis_dat %>% mutate(country=ifelse(country_region=="Mainland China", "China", country))
  arcgis_dat <- arcgis_dat %>% mutate(country=ifelse(province %in% c("Taiwan","Hong Kong","Tibet"), "China",country))
  arcgis_dat <- arcgis_dat %>% mutate(country=ifelse(is.na(country), as.character(province), country))
  arcgis_dat <- arcgis_dat %>% filter(country == "China") %>% select(-country_region)
  arcgis_dat <- arcgis_dat %>% mutate(n = ifelse(n < 0, 0, n))
  arcgis_dat <- expand_arcgis(arcgis_dat)
  
  subset_combined_dat <- combined_dat %>% filter(date_confirmation <= convert_date(switch_date) & country == "China")
  
  final <- bind_rows(arcgis_dat, subset_combined_dat) %>% arrange(province, date_confirmation)
  final
}

convert_datestring_to_date <- function(datestring){
  dates <- sapply(datestring, function(x){
    tmp <- strsplit(x, split="_")
    day <- as.numeric(str_extract(tmp[[1]][1], "[0-9]+"))
    month <- str_extract(tmp[[1]][1], "[aA-zZ]+")
    ## Extract time
    time <- str_extract(tmp[[1]][2], "[0-9]+")
    ## Only want first digit
    if(length(strsplit(time, "")[[1]]) == 3){
      time <- as.numeric(strsplit(time, "")[[1]][1])
    } else if(length(strsplit(time, "")[[1]]) == 4) {
      time <- as.numeric(paste(strsplit(time, "")[[1]][1:2], collapse= ""))
    } else {
      time <- as.numeric(time)
    }
    ## AM or PM?
    am <- str_extract(tmp[[1]][2], "[aA-zZ]+")
    if((am == "PM" | am == "pm") & time != 12) {
      time <- time + 12
    }
    if((am == "am" | am == "AM") & time == 12) {
      time <- 0
    }
    final_time <- make_datetime(2020, pmatch(month, month.name), day, time, 0)
    if(hour(final_time) == 0) final_time <- final_time - minutes(1)
    final_time
    },
    simplify=FALSE)
  do.call("c", dates)
}


## The updated ARC GIS dates mess with my pipeline, so for now manually plugging in the dates
sheet_name_key <- c("Feb02_9PM"="2020-02-02 21:00:00 UTC", 
                    "Feb02_745pm"="2020-02-02 19:45:00 UTC", 
                    "Feb02_5am"="2020-02-02 05:00:00 UTC",  
                    "Feb01_11pm"="2020-02-01 23:00:00 UTC",  
                    "Feb01_6pm"="2020-02-01 18:00:00 UTC",  
                    "Feb01_10am"="2020-02-01 10:00:00 UTC",  
                    "Jan31_7pm"="2020-01-31 19:00:00 UTC",  
                    "Jan31_2pm"="2020-01-31 14:00:00 UTC",  
                    "Jan30_930pm"="2020-01-30 21:30:00 UTC",  
                    "Jan30_11am"="2020-01-30 11:00:00 UTC",  
                    "Jan29_9pm"="2020-01-29 21:00:00 UTC",  
                    "Jan29_230pm"="2020-01-29 14:30:00 UTC",  
                    "Jan29_130pm"="2020-01-29 13:30:00 UTC",  
                    "Jan28_11pm"="2020-01-28 23:00:00 UTC",  
                    "Jan28_6pm"="2020-01-28 18:00:00 UTC",  
                    "Jan28_1pm"="2020-01-28 13:00:00 UTC", 
                    "Jan27_830pm"="2020-01-27 20:30:00 UTC", 
                    "Jan27_7pm"="2020-01-27 19:00:00 UTC",  
                    "Jan27_9am"="2020-01-27 09:00:00 UTC",  
                    "Jan26_11pm"="2020-01-26 23:00:00 UTC",  
                    "Jan26_11am"="2020-01-26 11:00:00 UTC", 
                    "Jan25_10pm"="2020-01-25 22:00:00 UTC",  
                    "Jan25_12pm"="2020-01-25 12:00:00 UTC",  
                    "Jan25_12am"="2020-01-24 23:59:59 UTC", 
                    "Jan24_12pm"="2020-01-24 12:00:00 UTC",  
                    "Jan24_12am"="2020-01-23 23:59:59 UTC",  
                    "Jan23_12pm"="2020-01-23 12:00:00 UTC",  
                    "Jan22_12pm"="2020-01-22 12:00:00 UTC",  
                    "Jan22_12am"="2020-01-21 23:59:59 UTC")
#' fit geometric distribution using Stan
fit_geometric_stan <- function(delay_data, model) {
  
  # define data and constants
  data_list <- list(  
    N= length(delay_data),    # number of sampling times
    delay = delay_data)
  
  # fit the model to data
  chains <- 3 # number of NUTS chains to run in parallel
  set.seed(1) # set random number generator seed for reproducibility
  # fit model to data using No U-Turn sampler
  if(missing(model)) {
    fit <- stan(get_model_filename("geometric.stan"), data=data_list,  chains=3,verbose = TRUE)
  } else {
    fit <- sampling(model, data=data_list,  chains=3,verbose = TRUE)
  }
  fit
}

#' simulate confirmation times
simulate_confirmation_times <- function(date_onset_symptoms, p_confirm_delay){
  # start from symptom onset time for now
  
  ## From the provided geometric distribution + 1
  sim_confirmation_delays <- rgeom(length(date_onset_symptoms), p_confirm_delay) + 1
  date_confirmation <- date_onset_symptoms + floor(sim_confirmation_delays)

  return(date_confirmation)
}

plot_forward_simulation <- function(data_quantiles, onset_data, max_date="27.01.2020",
                                min_date="01.12.2019", cols = c("grey40"),
                                cols2 = c("grey40"),
                                ymax=500,ybreaks=25, 
                                title = "Forward simulation of confirmed cases from symptom onsets"){
  
  p <- ggplot(data_quantiles)
  p <- p +
    #geom_rect(data=threshold_dat,aes(xmin=xmin,xmax=xmax,ymin=0,ymax=ymax,alpha=fills),fill="red") +
    geom_bar(data=onset_data,aes(x=date_onset_symptoms,y=n),fill = "orange", stat="identity") +
    geom_ribbon(aes(x=date,ymax=upper,ymin=lower),alpha=0.25, fill = "grey40") +
    geom_line(aes(x=date, y=mean),size=1, col = "grey40") +
    scale_y_continuous(expand=c(0,0),breaks=seq(0,ymax,by=ybreaks)) +
    coord_cartesian(ylim=c(0,ymax),xlim=c(convert_date(min_date), convert_date(max_date)+1)) +
    scale_x_date(limits=c(convert_date("01.12.2019"),convert_date(max_date)+1),
                 breaks="5 day") + 
    scale_fill_manual(values= cols) + 
    scale_color_manual(values=cols2,guide="none") +
    ggtitle(title) +
    ylab("Count") + xlab("Date of event") +
    theme_pubr() +
    theme(axis.text.x=element_text(angle=45,hjust=1),
          panel.grid.major = element_line(colour="grey70"),
          legend.position = c(0.25,0.75)) 
  p
}
>>>>>>> bc58d52dc70cbea06ed0ab0bd01fba1aa1db27a7
