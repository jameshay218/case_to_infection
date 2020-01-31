## Get median start date for hubei
min_date_overall <- province_data %>% filter(n > 0) %>% 
  filter(province=="Hubei") %>% group_by(repeat_no) %>% 
  summarise(start_hubei=min(date)) %>% ungroup() %>% summarise(y=median(start_hubei)) %>% pull(y)

## Get median start dates for all provinces
## Get how much you need to shift this to give the 
## same value as for Hubei
min_dates <- province_data %>% filter(n > 0) %>% group_by(province,repeat_no) %>% 
  summarise(min_prov=min(date)) %>% ungroup() %>% group_by(province) %>% summarise(start_prov=median(min_prov))
min_dates$hubei_start <- min_date_overall

min_dates <- min_dates %>% 
  mutate(shift_needed = as.numeric(hubei_start-start_prov))

## Merge these and then generate a variable that shifts 
## confirmation days to give the same median
## as Hubei
province_data1 <- province_data %>% left_join(min_dates)
province_data1 <- province_data1 %>% filter(n > 0) 
province_data1 <- province_data1 %>% filter(var=="date_infection") %>% mutate(x=as.numeric(date +shift_needed))
#province_data1$x <- convert_date(province_data1$x-province_data1$hubei_start) 

## Now get quantiles on this
#province_data1 <- province_data1 %>% group_by(province, repeat_no) %>% 
#  summarise(start=min(x)) %>% 
#  right_join(province_data1) %>% mutate(x=x-start)
province_data1$n_inflated_log <- log10(province_data1$n_inflated)
sim_data_quantiles_province <- province_data1 %>% ungroup() %>% group_by(x, province) %>% 
  do(data.frame(t(c(quantile(.$n_inflated_log, probs = c(0.025,0.5,0.975),na.rm=TRUE),mean=mean(.$n_inflated))))) %>% ungroup()


colnames(sim_data_quantiles_province) <- c("date","province","lower","median","upper","mean")
factor_order <- as.character(total_confirmed_prov$province)

sim_data_quantiles_province$province <- factor(as.character(sim_data_quantiles_province$province), 
                                               levels=factor_order)
subset_hubei <- sim_data_quantiles_province %>% filter(province == "Hubei") %>% select(date, lower, median, upper)
colnames(subset_hubei) <- c("date","lower_h","median_h","upper_h")
dat <- left_join(sim_data_quantiles_province, subset_hubei)
#max_help <- dat %>%
#  group_by(province) %>%
#  summarise(max_truncate=max(median)+5) %>% right_join(dat) %>%
#  mutate(lower_h=ifelse(lower_h >= max_truncate,  max_truncate, lower_h)) %>%
#  mutate(upper_h=ifelse(upper_h >=  max_truncate,  max_truncate, upper_h)) %>%
#  mutate(median_h=ifelse(median_h >  max_truncate,  NA, median_h)) %>% 
#  ungroup()
max_help <- dat  
max_help$date <- max_help$date - min(max_help$date)
p_comparison <- ggplot(max_help) + 
  geom_line(aes(x=date,y=median_h),col=orange_color,size=0.75) + 
  geom_ribbon(aes(x=date, ymin=lower_h,ymax=upper_h),alpha=0.5,fill=orange_color) +
  geom_ribbon(aes(x=date,ymin=lower,ymax=upper),alpha=0.5,fill=blue_color) +
  geom_line(aes(x=date,y=median),col=blue_color,size=0.75) + 
  coord_cartesian(xlim=c(80,170)) +
  ylab("log10 infection incidence") +
  xlab("Days since median epidemic start day for that province") +
  theme_bw() +
  #geom_vline(xintercept=0,linetype="dashed") +
  #geom_hline(yintercept=0, linetype="dashed",col="red") #+
  facet_wrap(~province, scales="free_y", ncol=5)
p_comparison

png("plots/comparison_plot.png",height=10,width=10,units="in",res=300)
plot(p_comparison)
dev.off()




