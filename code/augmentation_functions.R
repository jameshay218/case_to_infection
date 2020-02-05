generate_forward_probabilities_dist <- function(repeats, all_prob_parameters, tmax=100){
  ## Make total confirmation delay distribution
  #ret <- expand_grid(repeat_no=1:repeats,confirm_delay=0:tmax,symp_delay=0:tmax) %>% 
  #  mutate(total_delay=confirm_delay+symp_delay) %>% full_join(weibull_pars) %>%
  #  mutate(total_prob=forward_prob_total(symp_delay, confirm_delay, alpha, sigma, p_confirm_delay)) %>%
  #  group_by(total_delay, repeat_no) %>% summarise(total_prob_sum=sum(total_prob)) %>% ungroup() %>% 
  #  group_by(repeat_no) %>%
  #  mutate(cumu_prob_total=cumsum(total_prob_sum))
    
  delay_table <- expand_grid(repeat_no=1:repeats,confirm_delay=0:tmax)
  confirmation_delay_table <- all_delay_prob_parameters %>% select(repeat_no, gamma_scale, gamma_shape, date_confirmation) %>%
    right_join(delay_table) %>% mutate(confirm_prob=ddgamma(confirm_delay-1, scale=gamma_scale, shape=gamma_shape)) %>%
    arrange(repeat_no, date_confirmation, confirm_delay) %>%
    group_by(repeat_no,date_confirmation) %>% mutate(cumu_prob_confirm=cumsum(confirm_prob)) %>% ungroup()
  
  
  ## Make just final confirmation delay distribution
  confirmation_delay_geometric <- all_delay_prob_parameters %>% select(repeat_no, fixed_geom) %>%
    distinct() %>%
    right_join(delay_table) %>%
    mutate(confirm_prob=dgeom(confirm_delay-1,fixed_geom)) %>%
    group_by(repeat_no) %>% mutate(cumu_prob_confirm=cumsum(confirm_prob)) %>% ungroup()
  
  symptom_delay_table <- expand_grid(repeat_no=1:repeats,symp_delay=0:tmax)
  
  ## Make just symptom onset delay distribution
  symptom_delay <- all_delay_prob_parameters %>% select(repeat_no, alpha, sigma) %>% 
    distinct() %>%
    right_join(symptom_delay_table) %>%
    mutate(symp_prob=incubation_prob(symp_delay, alpha, sigma)) %>%
    group_by(repeat_no) %>% mutate(cumu_prob_symp=cumsum(symp_prob)) %>% ungroup()
  
  return(list(confirmation_delay_table, confirmation_delay_geometric, symptom_delay))
}

generate_total_forward_probabilities_dist <- function(repeats, weibull_pars, p_confirm_delay, tmax=100){
  ## Make total confirmation delay distribution
  ret <- expand_grid(repeat_no=1:repeats,confirm_delay=0:tmax,symp_delay=0:tmax) %>% 
    mutate(total_delay=confirm_delay+symp_delay) %>% full_join(weibull_pars) %>%
    mutate(total_prob=forward_prob_total(symp_delay, confirm_delay, alpha, sigma, p_confirm_delay)) %>%
    group_by(total_delay, repeat_no) %>% summarise(total_prob_sum=sum(total_prob)) %>% ungroup() %>% 
    group_by(repeat_no) %>%
    mutate(cumu_prob_total=cumsum(total_prob_sum))
  
  return(ret)
}


incubation_prob  <- function(symp_delay, alpha, sigma){
  pweibull(symp_delay+1, alpha, sigma) - pweibull(symp_delay, alpha, sigma)
}
forward_prob_total <- function(symp_delay, confirm_delay, alpha, sigma, p_confim_delay){
  (pweibull(symp_delay+1, alpha, sigma) - pweibull(symp_delay, alpha, sigma))*dgeom(confirm_delay-1, p_confirm_delay)
  
}



## Generates distributions for the probability of an infection/symptom onset occuring on a given day being
## observed at every day in the future
generate_forward_probabilities <- function(alpha, sigma, p_confirm_delay,tmax=100){
  ret <- expand_grid(confirm_delay=0:tmax,symp_delay=0:tmax) %>% 
    mutate(total_delay=confirm_delay+symp_delay) %>% 
    mutate(total_prob=forward_prob_total(symp_delay, confirm_delay, alpha, sigma, p_confirm_delay)) %>%
    group_by(total_delay) %>% summarise(total_prob_sum=sum(total_prob)) %>% ungroup() %>% 
    mutate(cumu_prob=cumsum(total_prob_sum))
  return(ret)
}


#' Augment infection times
augment_infection_times <- function(dat, inc_period_alpha, inc_period_sigma, p_confirm_delay=0.15, minimum_confirmation_delay=1){
  ## Sim symptom onset times for those without
  which_to_sim <- which(is.na(dat$date_onset_symptoms) & !is.na(dat$date_confirmation))
  n_to_sim <- nrow(dat[which_to_sim,])
  
  ## From the provided geometric distribution + 1
  if(!is.null(p_confirm_delay)) {
    sim_confirmation_delays <- rgeom(n_to_sim, p_confirm_delay) + minimum_confirmation_delay
  } else {
    sim_confirmation_delays <- dat %>% mutate(confirmation_delay=rdgamma(n(), scale=gamma_scale_backward, shape=gamma_shape_backward)+1) %>% 
      pull(confirmation_delay)
  }
  dat[which_to_sim,"date_onset_symptoms"] <- dat[which_to_sim,"date_confirmation"] - sim_confirmation_delays
  
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

plot_augmented_events <- function(data_quantiles, confirmed_data, 
                                  var_name="date_infections",
                                      max_date="01.02.2020",min_date="01.12.2019",
                                      ymax=500,ybreaks=25,
                                      thresholds=NULL,
                                      cols=c("skyblue","blue"),
                                      title="Augmented and observed infection incidence against confirmed cases",
                                      var_labels=c("95th percentile on observed and unobserved infections (estimates)",
                                               "Number of infections from observed cases (estimates)"),bot=TRUE){
  data_quantiles[data_quantiles$inflated == "total",c("median","mean")] <- NA
  data_quantiles$inflated <- factor(data_quantiles$inflated, levels=c("total","inflated0"))
  p <- ggplot(data_quantiles[data_quantiles$var == var_name,])
  
  ## Threshold markers for infection onsets
  if(!is.null(thresholds)) {
    threshold_dat <- data.frame(text=c(">95%",">80%",">50%",">20%"),
                                x_val=thresholds,stringsAsFactors = FALSE)
    p <- p +    
      geom_vline(xintercept=thresholds,linetype="dashed") +
      geom_label(data=threshold_dat,aes(x=x_val, y=ymax*0.8, label=text))
  }
  p <- 
    p +
    geom_bar(data=confirmed_data,aes(x=date_confirmation,y=n,fill=var),stat="identity",col="black", size=0.5) +
    geom_ribbon(aes(x=date,ymax=upper,ymin=lower,fill=inflated),alpha=0.4, size=0.5) +
    geom_line(aes(x=date, y=mean,col=inflated),size=0.5) +
    scale_y_continuous(expand=c(0,0),breaks=seq(0,ymax,by=ybreaks)) +
    coord_cartesian(ylim=c(0,ymax),xlim=c(convert_date(min_date), convert_date(max_date)+1)) +
    scale_x_date(limits=c(convert_date(min_date),convert_date(max_date)+1),
                 breaks="5 day") + 
    scale_fill_manual(values= c("inflated0"=cols[2],"total"=cols[1],"confirmed"="grey40"), name="Variable",
                      labels=c("Confirmed cases (reports to date)", var_labels[c(2,1)]))+ 
    scale_color_manual(values= c("inflated0"=cols[2],"total"=cols[1],"confirmed"="grey40"), guide="none") +
    ylab("Count") + xlab("Date of event") +
    theme_pubr()  +
    ggtitle(title)
  if(bot) {
    p <- p +
    theme(axis.text.x=element_text(angle=45,hjust=1),
          panel.grid.major = element_line(colour="grey70",size=0.2),
          legend.position = c(0.3,0.6))
  } else {
    p <- p + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank(),
                   panel.grid.major = element_line(colour="grey70",size=0.2),
                   legend.position = c(0.3,0.6))
  }
  p
}

plot_augmented_events_byprovince <- function(data_quantiles_province, confirmed_data_province, 
                                             provinces=NULL,
                                  var_name="date_infections",
                                  max_date="01.02.2020",min_date="01.12.2019",
                                  thresholds=NULL,
                                  ncol=5,
                                  cols=c("skyblue","blue"),
                                  title="Augmented and observed infection incidence against confirmed cases",
                                  var_labels=c("95th percentile on observed and unobserved infections (estimates)",
                                               "Number of infections from observed cases (estimates)")){
  data_quantiles_province[data_quantiles_province$inflated == "total",c("median","mean")] <- NA
  data_quantiles_province$inflated <- factor(data_quantiles_province$inflated, levels=c("total","inflated0"))
  data_quantiles_province <- data_quantiles_province %>% filter(province %in% provinces)
  confirmed_data_province <- confirmed_data_province %>% filter(province %in% provinces)
  
  p <- ggplot(data_quantiles_province[data_quantiles_province$var == var_name,])
  p <- p +
    geom_bar(data=confirmed_data_province,aes(x=date_confirmation,y=n,fill=var),stat="identity",col="black", size=0.5) +
    geom_ribbon(aes(x=date,ymax=upper,ymin=lower,fill=inflated),alpha=0.4, size=0.5) +
    geom_line(aes(x=date, y=median,col=inflated),size=0.5) +
    scale_y_continuous(expand=c(0,0)) +
    coord_cartesian(xlim=c(convert_date(min_date), convert_date(max_date)+1)) +
    scale_x_date(limits=c(convert_date(min_date),convert_date(max_date)+1),
                 breaks="7 day") + 
    facet_wrap(~province, scales="free_y",ncol=ncol) +
    scale_fill_manual(values= c("inflated0"=cols[2],"total"=cols[1],"confirmed"="grey40"), name="Variable",
                      labels=c("Confirmed cases (reports to date)", var_labels[c(2,1)]))+ 
    scale_color_manual(values= c("inflated0"=cols[2],"total"=cols[1],"confirmed"="grey40"), guide="none") +
    ylab("Count") + xlab("Date of event") +
    theme_pubr()  +
    ggtitle(title) +
      theme(axis.text.x=element_text(angle=45,hjust=1),
            panel.grid.major = element_line(colour="grey70",size=0.2),
            legend.position ="none")
  p
}



plot_augmented_data <- function(data_quantiles, confirmed_data, 
                                max_date="27.01.2020",
                                min_date="01.12.2019", 
                                cols1 = c("skyblue","blue"),
                                cols2 = c("orange","red"),
                                ymax1=500,
                                ymax2=500,
                                ybreaks=25, 
                                thresholds=NULL,
                                thresholds_symp=NULL,
                                title1 = "Augmented and observed timings of infection incidence in China",
                                title2 = "Augmented and observed timings of symptom onset incidence in China",
                                var_labels1=c("95th percentile on observed and yet-to-be\n observed infections (estimates)",
                                             "Number of infections from observed cases (estimates)"),
                                var_labels2=c("95th percentile on observed and yet-to-be\n observed symptom onsets (estimates)",
                                             "Number of symptom onsets from observed cases (estimates)")){

  p1 <- plot_augmented_events(data_quantiles, confirmed_data, var_name="date_infections",
                              max_date, min_date, ymax1, ybreaks, thresholds, cols1, title1,
                              var_labels1,bot=FALSE)
  p2 <- plot_augmented_events(data_quantiles, confirmed_data, var_name="date_onset_symptoms",
                              max_date, min_date, ymax2, ybreaks, thresholds_symp, cols2, title2,
                              var_labels2,bot=TRUE)
  #plot_grid(p1,p2,ncol=1,align="h")
  p1 / p2
}

expand_arcgis <- function(arcgis_dat) {
  arcgis_dat %>% group_by(date_confirmation, province, country) %>% filter(n > 0) %>% drop_na() %>% uncount(n) %>% ungroup()
}

merge_data <- function(linelist_dat, arcgis_dat, switch_date){
  arcgis_dat <- arcgis_dat %>% select(province, raw_day, country_region, country, diff)
  colnames(arcgis_dat) <- c("province","date_confirmation","country_region","country","n")
  arcgis_dat <- arcgis_dat %>% mutate(country=ifelse(country_region=="Mainland China", "China", country))
  arcgis_dat <- arcgis_dat %>% mutate(country=ifelse(province %in% c("Taiwan","Hong Kong","Tibet"), "China",country))
  arcgis_dat <- arcgis_dat %>% mutate(country=ifelse(is.na(country), as.character(province), country))
  arcgis_dat <- arcgis_dat %>% filter(country == "China") %>% select(-country_region)
  arcgis_dat <- arcgis_dat %>% mutate(n = ifelse(n < 0, 0, n))
  arcgis_dat <- expand_arcgis(arcgis_dat) %>% select(province,date_confirmation) %>% 
    mutate(province = as.character(province))
  arcgis_dat <- arcgis_dat %>% mutate(date_onset_symptoms=NA, 
                                      date_death_or_discharge=NA,
                                      date_admission_hospital=NA)
  
  
  subset_combined_dat <- linelist_dat %>% 
    filter(date_confirmation <= convert_date(switch_date) & country == "China") %>%
    select(province, date_confirmation,date_death_or_discharge,date_onset_symptoms,date_admission_hospital) %>%
    mutate(province=as.character(province))
  
  final <- bind_rows(arcgis_dat, subset_combined_dat) %>% arrange(province, date_confirmation) %>% ungroup()
  final
}


#' simulate confirmation times
simulate_confirmation_times <- function(date_onset_symptoms, p_confirm_delay=0.15, gamma_mean=NULL, gamma_var=NULL,
                                        minimum_confirmation_delay=1){
  # start from symptom onset time for now
  ## From the provided geometric distribution + 1
  if (!is.null(gamma_mean) & !is.null(gamma_var)){
    sim_confirmation_delays <- ceiling(rgamma_mean(length(date_onset_symptoms), gamma_mean, gamma_var))
  } else {
    sim_confirmation_delays <- rgeom(length(date_onset_symptoms), p_confirm_delay) + minimum_confirmation_delay
  }
  date_confirmation <- date_onset_symptoms + sim_confirmation_delays
  
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
          panel.grid.major = element_line(colour="grey70",size=0.2),
          legend.position = c(0.25,0.75)) 
  p
}
