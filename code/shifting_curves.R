## Get median start date for hubei
min_date_overall <- province_data %>% filter(n > 0) %>% 
  filter(province=="Hubei") %>% group_by(repeat_no) %>% 
  summarise(x=min(date)) %>% ungroup() %>% summarise(y=median(x)) %>% pull(y)

## Get median start dates for all provinces
## Get how much you need to shift this to give the 
## same value as for Hubei
min_dates <- province_data %>% filter(n > 0) %>% group_by(province,repeat_no) %>% 
  summarise(min_prov=min(date)) %>% ungroup() %>% group_by(province) %>%
  summarise(min_prov_med=median(min_prov)) %>%
  mutate(offset=min_prov_med-min_date_overall) 

## Merge these and then generate a variable that shifts 
## confirmation days to give the same median
## as Hubei
province_data1 <- province_data %>% left_join(min_dates)
province_data1 <- province_data1 %>% filter(n > 0) 
province_data1 <- province_data1 %>% filter(var=="date_infection") %>% mutate(x=as.numeric(date - offset))
province_data1$x <- convert_date(province_data1$x) 
## Start at 0
## province_data1$x1 <- province_data1$x - min(province_data1$x)

wow <- province_data1 %>% ungroup() %>% 
  select(repeat_no, var, x, province, n) %>%
  group_by(province) %>% 
  complete(repeat_no, var, x, province, fill=list(n=0))

province_data1$x1 <- as.numeric(province_data1$min_prov_med - province_data1$x)
## Find Hubei's incidence for each day for each repeat
ref_points <- province_data1 %>% filter(province=="Hubei" & var == "date_infection") %>% select(repeat_no, x1, n)
colnames(ref_points)[3] <- "hubei_n"

## Merge these references with the data per province
wowow <- province_data1 %>% filter(var == "date_infection") %>% select(repeat_no, x1, province, n) %>% left_join(ref_points)
## And find difference between hubei and each other location on that day
wowow$diffn <- wowow$n - wowow$hubei_n
wowow %>% filter(province=="Hubei" & n > 0) %>% View


## Now get quantiles on this
sim_data_quantiles_province <- wowow %>% ungroup() %>% group_by(x1, province) %>% 
  do(data.frame(t(quantile(.$diffn, probs = c(0.025,0.5,0.975),na.rm=TRUE))))


colnames(sim_data_quantiles_province) <- c("date","province","lower","median","upper")
factor_order <- as.character(total_confirmed_prov$province)

sim_data_quantiles_province$province <- factor(as.character(sim_data_quantiles_province$province), 
                                               levels=factor_order)

ggplot(sim_data_quantiles_province) + 
  geom_ribbon(aes(x=date,ymin=lower,ymax=upper),alpha=0.25) +
  geom_line(aes(x=date,y=median)) + 
  ylab("Cases in province on day X minus cases in Hubei on day X") +
  xlab("Day X (days since median epidemic start day)") +
  theme_bw() +
  geom_vline(xintercept=0,linetype="dashed") +
  geom_hline(yintercept=0, linetype="dashed",col="red") +
  facet_wrap(~province, scales="free_y")





