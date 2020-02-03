## Pull latest cumulative incidence data from arcgis
url <- "https://docs.google.com/spreadsheets/d/1UF2pSkFTURko2OvfHWWlFpDFAr1UxCBA4JLwlSP6KFo/htmlview?usp=sharing&sle=true"

colnames_key <- c("Province/State"="province", "Country/Region"="country_region", "Last Update"="updated", 
                 "Confirmed"="confirmed","Deaths"="deaths", "Recovered"="recovered","Demised"="deaths",
                 "Date last updated"="updated","Suspected"="suspected",
                 "Country"="country_region")

all_dat <- read_sheet(ss=url)
all_dat <- all_dat %>% mutate(`Province/State` = ifelse(is.na(`Province/State`), `Country/Region`,`Province/State`))
all_dat <- all_dat %>% pivot_longer(-c("Province/State","Country/Region","First confirmed date in country (Est.)","Lat","Long"),
                                    names_to="updated",values_to="confirmed") 
colnames(all_dat) <- c("province","country_region","first_confirmed","lat","long","updated","confirmed")
use_data <- all_dat
use_data$updated_day <- sapply(use_data$updated, function(x) strsplit(x, split=" ")[[1]][[1]])
use_data$updated_day <- convert_date(use_data$updated_day,format="%m/%d/%Y")
use_data$updated <- strptime(use_data$updated, format="%m/%d/%Y %I:%M %p")
use_data$updated <- ymd_hms(use_data$updated)
## The other factors are in days, so can only use daily incidence.
## I reason that we should use the last report of each day that exists for
## each province
## Get the last report for each day
subset_report_dates <-  use_data %>% select(updated, updated_day, province) %>%  
  arrange(province,updated_day, updated) %>% 
  group_by(updated_day, province) %>%
  unique %>% mutate(n=1:n()) %>% 
  group_by(province, updated_day) %>%
  filter(n==max(n)) %>%
  arrange(province, updated_day) %>% 
  ungroup()


## Get date as day rather than POSIT
use_data_subset <- use_data %>% right_join(subset_report_dates)
use_data_subset <- use_data_subset %>% mutate(raw_day=as.Date(updated,origin="1970-01-01")) 
use_data_subset <- use_data_subset %>% complete(province, raw_day, fill=list(confirmed=0))
use_data <- use_data %>% mutate(raw_day=as.Date(updated,origin="1970-01-01"))

## Reorder provinces so highest inc is first
use_data_subset$province <- fct_reorder(use_data_subset$province, use_data_subset$confirmed,.fun=max,.desc = TRUE)
use_data$province <- fct_reorder(use_data$province, use_data$confirmed,.fun=max,.desc = TRUE)

p_cumu_confirmed_perday <- ggplot(use_data_subset) +
  geom_bar(aes(x=raw_day,y=confirmed),stat="identity") + 
  facet_wrap(~province,scales="free_y",ncol=5)+
  theme_bw() +
  scale_x_date(limits=c(convert_date("19.01.2020"),convert_date("31.01.2020")),
               breaks="1 day") +
  xlab("Day") + ylab("Cumulative confirmed cases per day") +
  theme(axis.text.x=element_text(angle=45,hjust=1))

png("plots/reports_cumu_conf_perday.png",width=12,height=12,res=300,units="in")
p_cumu_confirmed_perday
dev.off()



p_cumu_confirmed_per_report <- ggplot(use_data[use_data$province == "Hubei",]) +
  geom_point(aes(x=updated,y=confirmed)) + 
  geom_point(data=use_data_subset[use_data_subset$province=="Hubei",],aes(x=as.POSIXct(updated_day),y=confirmed),col="red") +
  facet_wrap(~province,scales="free_y",ncol=5) +
  scale_x_datetime(date_breaks="1 day") +
  theme_bw() +
  xlab("Time of report") + ylab("Cumulative confirmed cases per time of report") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) 

png("plots/reports_cumu_conf_perreport.png",width=12,height=12,res=300,units="in")
p_cumu_confirmed_per_report
dev.off()



use_data_diff <- use_data_subset %>% 
  group_by(province) %>%
  mutate(prev=dplyr::lag(confirmed,n=1)) %>%
  mutate(prev=ifelse(raw_day==min(raw_day), NA, prev)) %>%
  mutate(diff=confirmed-prev) 

use_data_diff$country <- use_data_diff$country_region

p_diff_per_day <- ggplot(use_data_diff) + 
    geom_bar(aes(x=raw_day,y=diff),stat="identity") + 
    facet_wrap(~province, scale="free_y",ncol=5) +
    theme_bw() +
    scale_x_date(limits=c(convert_date("19.01.2020"),convert_date("31.01.2020")),
                 breaks="1 day") +
    xlab("Day of report") + ylab("New confirmed cases from previous day") +
    theme(axis.text.x=element_text(angle=45,hjust=1)) 
png("plots/reports_new_perday.png",width=12,height=12,res=300,units="in")
p_diff_per_day
dev.off()


