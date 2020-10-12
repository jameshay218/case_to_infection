#####
## NOTE - NEED TO RUN SIMPLER_METHOD.R FIRST

## Plot example province
date_min <- as.Date("2020-01-01")
date_max <- as.Date('2020-03-02')
use_prov <- "Guangdong"
use_city <- "Guangzhou"
pA <- all_incidence_province %>% 
  filter(province_raw == use_prov,
         dates >= date_min & dates <= date_max) %>%
  ggplot() +
  geom_bar(data=confirmed_cases_date %>% 
             filter(province_raw==use_prov,
                    dates >= date_min & dates <= date_max),
           aes(x=dates,y=n),stat="identity",col="grey5", size=0.5,fill="grey50") +
  geom_ribbon(aes(x=dates,ymax=n_onset,ymin=0),alpha=0.2,fill="#56B4E9") +
  geom_line(aes(x=dates,y=n_onset),col="#56B4E9",size=0.75) +
  geom_ribbon(aes(x=dates,ymax=n_infected,ymin=0),alpha=0.2,fill="#E69F00") +
  geom_line(aes(x=dates,y=n_infected),col="#E69F00",size=0.75) +
  scale_y_continuous(expand=c(0,0),limits=c(0,150),breaks=seq(0,150,by=25)) +
  scale_x_date(limits=c(date_min, date_max+5),labels=seq(date_min, date_max+5,by="7 days"),breaks=seq(date_min, date_max+5,by="7 days")) +
  xlab("Date") +
  ylab("Incidence (number of cases)") +
  export_theme+
  ggtitle("Guangdong") +
  theme(#axis.text.x=element_blank(),
        #axis.ticks.x = element_blank(),
        axis.title.x=element_blank(),
        plot.margin = margin(0,0,0,0, "cm")) + 
  labs(tag="A")

pB <- prov_inc_prev_cali %>% filter(city == use_city) %>%
  ggplot() + 
  geom_line(aes(x=date,y=travel_prev),col="grey5",size=0.75) + 
  geom_ribbon(aes(x=date,ymax=travel_prev,ymin=0),fill="#009E73") + 
  geom_line(data=city_n_inf_caladj_den %>% filter(city == use_city), 
            aes(x=date,y=n_infected_caladj),col="grey5",size=0.75) +
  geom_ribbon(data=city_n_inf_caladj_den %>% filter(city == use_city),
              aes(x=date,ymax=n_infected_caladj,ymin=0),fill="#E69F00") +
  scale_x_date(limits=c(date_min, date_max+5),labels=seq(date_min, date_max+5,by="7 days"),breaks=seq(date_min, date_max+5,by="7 days")) +
  scale_y_continuous(expand=c(0,0),limits=c(0,0.00002),breaks=seq(0,0.00002,by=0.000005)) +
  ylab("Prevalence/incidence (per capita)") +
  xlab("Date") +
  ggtitle("Guangzhou") +
  export_theme +
  theme(#axis.text.x=element_blank(),
        #axis.ticks.x = element_blank(),
        axis.title.x=element_blank(),
        plot.margin = margin(0,0,0,0, "cm")) + 
  labs(tag="B")

pA | pB

#mt <- read_csv("./data/master_table.csv",guess_max = Inf)
prev_all <- mt %>% select(is_wuhan, prevalence_o, date, scenario,origin_city) %>% distinct()

prev_summary <- prev_all %>% 
  group_by(date, is_wuhan, scenario) %>% 
  summarize(prev_average=mean(prevalence_o)) %>%
  rename(Scenario=scenario)

pC <- prev_summary %>% filter(is_wuhan ==1) %>%
  ggplot() + 
  geom_line(aes(x=date,y=prev_average,col=Scenario)) +
  #geom_jitter(aes(x=date,y=prev_average,col=Scenario)) +
  scale_x_date(limits=c(date_min, date_max+5),labels=seq(date_min, date_max+5,by="7 days"),breaks=seq(date_min, date_max+5,by="7 days")) +
  scale_y_continuous(expand=c(0,0),limits=c(0,0.005),breaks=seq(0,0.005,by=0.001)) +
  export_theme +
  ylab("Prevalence (per capita)") +
  xlab("Date") +
  theme(legend.position=c(0.8,0.7),
        legend.title=element_blank(),
        plot.margin = margin(0,0,0,0, "cm"))+
  ggtitle("Wuhan") +
  labs(tag="C")
pD <- prev_summary %>% filter(is_wuhan ==0) %>%
  ggplot() + geom_line(aes(x=date,y=prev_average,col=Scenario))+
  scale_x_date(limits=c(date_min, date_max+5),labels=seq(date_min, date_max+5,by="7 days"),breaks=seq(date_min, date_max+5,by="7 days")) +
  scale_y_continuous(expand=c(0,0),limits=c(0,0.0005),breaks=seq(0,0.0005,by=0.0001)) +
  xlab("Date") +
  ylab("Prevalence (per capita)") +
  ggtitle("China (excluding Wuhan)") +
  export_theme +
  theme(legend.position=c(0.8,0.7),
        legend.title=element_blank(),
        plot.margin = margin(0,0,0,0, "cm"))+
  labs(tag="D")


pdf("figures/Fig1.pdf",height=5,width=8)
(pA | pB)/ (pC | pD)
dev.off()
