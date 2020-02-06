## How many cases needed to be considered "started"?
cases_needed <- 1
## Get median start date for hubei
min_date_overall <- final_all_province %>% filter(total > cases_needed) %>% 
  filter(province=="Hubei") %>% group_by(repeat_no) %>% 
  summarise(start_hubei=min(date)) %>% ungroup() %>% summarise(y=median(start_hubei)) %>% pull(y)

## Get median start dates for all provinces
## Get how much you need to shift this to give the 
## same value as for Hubei
min_dates <- final_all_province %>% filter(total > cases_needed) %>% group_by(province,repeat_no) %>% 
  summarise(min_prov=min(date)) %>% ungroup() %>% group_by(province) %>% summarise(start_prov=median(min_prov))
min_dates$hubei_start <- min_date_overall

min_dates <- min_dates %>% 
  mutate(shift_needed = as.numeric(hubei_start-start_prov))

## Merge these and then generate a variable that shifts 
## confirmation days to give the same median start day as Hubei
## Get the full augmented data by province
province_data1 <- final_all_province %>% left_join(min_dates)
## Only interested in dates where more than threshold cases
province_data1 <- province_data1 %>% filter(total > cases_needed) 
## IMPORTANT - choose the threshold to cut off measurements from
print(threshold_vals)
print(thresholds)
threshold_i <- 1
province_data1 <- province_data1 %>% filter(date <= thresholds[threshold_i])
province_data1 <- province_data1 %>% filter(var=="date_infections") %>% mutate(x=as.numeric(date +shift_needed))


## Interested in log10 scale comparison
province_data1$n_inflated_log <- log10(province_data1$total)
sim_data_quantiles_province <- province_data1 %>% ungroup() %>% group_by(x, province) %>% 
  do(data.frame(t(c(quantile(.$n_inflated_log, probs = c(0.025,0.5,0.975),na.rm=TRUE),mean=mean(.$n_inflated_log))))) %>% 
  ungroup()


colnames(sim_data_quantiles_province) <- c("date","province","lower","median","upper","mean")
factor_order <- as.character(total_confirmed_prov$province)

sim_data_quantiles_province$province <- factor(as.character(sim_data_quantiles_province$province), 
                                               levels=factor_order)
subset_hubei <- sim_data_quantiles_province %>% filter(province == "Hubei") %>% select(date, lower, median, upper)
colnames(subset_hubei) <- c("date","lower_h","median_h","upper_h")
sim_data_quantiles_province$date <- sim_data_quantiles_province$date  - min(subset_hubei$date)
subset_hubei$date <- subset_hubei$date - min(subset_hubei$date)

sim_data_quantiles_province1 <- sim_data_quantiles_province %>% filter(province %in% factor_order[1:15])

p_comparison <- ggplot(sim_data_quantiles_province1) + 
  geom_line(data=subset_hubei, aes(x=date,y=median_h, col="#D55E00"),size=0.75) + 
  geom_ribbon(data=subset_hubei, aes(x=date, ymin=lower_h,ymax=upper_h, fill="#D55E00"),alpha=0.5) +
  geom_ribbon(aes(x=date,ymin=lower,ymax=upper,fill="#0072B2"),alpha=0.5) +
  geom_line(aes(x=date,y=median,col="#0072B2"),size=0.75) + 
  coord_cartesian(xlim=c(0,75)) +
  ylab("log10 infection incidence") +
  xlab("Days since median epidemic start day for that province") +
  ggtitle(paste0("Comparison of infection incidence across provinces as of ", thresholds[threshold_i])) +
  scale_fill_identity(name="",guide="legend",labels=c("Province","Hubei")) +
  scale_color_identity(name="",guide="legend",labels=c("Province","Hubei")) +
  theme_pubr()  +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        panel.grid.major = element_line(colour="grey70",size=0.2),
        legend.position =c(0.9,0.1),
        axis.title = element_text(size=12)) +
  facet_wrap(~province, scales="free_y", ncol=4)
p_comparison




