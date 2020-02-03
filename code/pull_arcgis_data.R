## Pull latest cumulative incidence data from arcgis
url <- "https://docs.google.com/spreadsheets/d/1wQVypefm946ch4XDp37uZ-wartW4V7ILdg-qYiDXUHM/htmlview?usp=sharing&sle=true#"
available_dates <- sheets_sheets(url)
#use_colnames <- c("Province/State", "Country/Region", "Last Update", "Confirmed", 
#                  "Deaths", "Recovered","Demised")

colnames_key <- c("Province/State"="province", "Country/Region"="country_region", "Last Update"="updated", 
                  "Last Update (UTC)"="updated",
                 "Confirmed"="confirmed","Deaths"="deaths", "Recovered"="recovered","Demised"="deaths",
                 "Date last updated"="updated","Suspected"="suspected",
                 "Country"="country_region")
use_colnames <- c("province","country_region","updated","confirmed","deaths","recovered")

## Go through and pull each sheet
all_data <- NULL
all_data_list <- NULL
for(date in available_dates) {
  Sys.sleep(0.01)
  tmp_dat <- read_sheet(ss=url,sheet=date,range="A:F")
  print(colnames(tmp_dat))
  tmp_dat <- tmp_dat[,colnames(tmp_dat) %in% names(colnames_key)]
  colnames(tmp_dat) <- colnames_key[colnames(tmp_dat)]
  use_colnames_tmp <- intersect(use_colnames, colnames(tmp_dat))
  tmp_dat <- tmp_dat %>% select(use_colnames_tmp)
  ## Some of the datetimes were entered incorrectly, so correct here
  if(date == "Jan25_12am") {
    tmp_dat <- tmp_dat %>% mutate(updated=as.POSIXct("2020-01-24 23:59:59 UTC",format="%Y-%m-%d %H:%M:%OS"))
  }
  if(date == "Jan25_10pm") {
    tmp_dat <- tmp_dat %>% mutate(updated=as.POSIXct("2020-01-25 22:00:00 UTC",format="%Y-%m-%d %H:%M:%OS"))
  } 
  if(date == "Jan24_12am") {
    tmp_dat <- tmp_dat %>% mutate(updated=as.POSIXct("2020-01-23 23:59:59 UTC",format="%Y-%m-%d %H:%M:%OS"))
  }
  colnames(tmp_dat)[colnames(tmp_dat) == "Demised"] <- "deaths"
  colnames(tmp_dat)[colnames(tmp_dat) == "Last Update (UTC)"] <- "updated"
  print(colnames(tmp_dat))
  tmp_dat$report_date <- date
  all_data <- bind_rows(all_data, tmp_dat)
  all_data_list[[date]] <- tmp_dat
}
use_data <- all_data
use_data$report_date <- convert_datestring_to_date(use_data$report_date)

use_data <- use_data %>% arrange(province, report_date)
use_data <- use_data %>% select(province, country_region,report_date,confirmed)
use_data <- unique(use_data)

## can't use data with no report date
use_data <- use_data %>% filter(!is.na(report_date))
## Replace province with smallest admin unit with entry
use_data <- use_data %>% mutate(province=ifelse(is.na(province), country_region, province)) %>%
  mutate(province=ifelse(is.na(province), country, province))
use_data$updated_day <- as.Date(use_data$report_date)
## The other factors are in days, so can only use daily incidence.
## I reason that we should use the last report of each day that exists for
## each province
## Get the last report for each day (sorry this is such a rubbish solution)
subset_report_dates <-  use_data %>% select(report_date, updated_day, province) %>%  
  arrange(province,updated_day, report_date) %>% 
  group_by(updated_day, province) %>%
  unique %>% mutate(n=1:n()) %>% 
  group_by(province, updated_day) %>%
  filter(n==max(n)) %>%
  arrange(province, updated_day) %>% 
  ungroup()

## Get date as day rather than POSIT
use_data_subset <- use_data %>% right_join(subset_report_dates)

## If contradictory results for the same update time, use lower cases
use_data_subset <- use_data_subset %>% 
  group_by(province, updated_day) %>% 
  arrange(province, updated_day, confirmed) %>% 
  mutate(n=1:n()) %>%
  group_by(province, updated_day) %>%
  filter(n==min(n)) %>%
  arrange(province, updated_day) %>% 
  ungroup()


use_data_subset <- use_data_subset %>% mutate(raw_day=as.Date(report_date,origin="1970-01-01")) 
use_data_subset <- use_data_subset %>% complete(province, raw_day, fill=list(confirmed=0))
use_data_subset <- use_data_subset %>% filter(!is.na(updated_day))
use_data <- use_data %>% mutate(raw_day=as.Date(report_date,origin="1970-01-01"))

## Reorder provinces so highest inc is first
use_data_subset$province <- fct_reorder(use_data_subset$province, use_data_subset$confirmed,.fun=max,.desc = TRUE)
use_data$province <- fct_reorder(use_data$province, use_data$confirmed,.fun=max,.desc = TRUE)

p_cumu_confirmed_perday <- ggplot(use_data_subset) +
  geom_bar(aes(x=raw_day,y=confirmed),stat="identity") + 
  facet_wrap(~province,scales="free_y",ncol=5)+
  theme_bw() +
  scale_x_date(limits=c(convert_date("19.01.2020"),convert_date(date_today)),
               breaks="1 day") +
  xlab("Day") + ylab("Cumulative confirmed cases per day") +
  theme(axis.text.x=element_text(angle=45,hjust=1))

png("plots/reports_cumu_conf_perday.png",width=12,height=12,res=300,units="in")
p_cumu_confirmed_perday
dev.off()



p_cumu_confirmed_per_report <- ggplot(use_data) +
  geom_point(aes(x=as.Date(report_date,origin="1970-01-01"),y=confirmed)) + 
  facet_wrap(~province,scales="free_y",ncol=5) +
  theme_bw() +
  scale_x_date(limits=c(convert_date("19.01.2020"),convert_date(date_today)),
             breaks="1 day") +
  xlab("Time of report") + ylab("Cumulative confirmed cases per time of report") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) 

png("plots/reports_cumu_conf_perreport.png",width=12,height=12,res=300,units="in")
p_cumu_confirmed_per_report
dev.off()


## Comparison of numbers from reports as of times reported
## and numbers used (last report of the day)
ggplot(use_data[use_data$province %in% c("Hubei","Guangdong","Zheijiang","Henan","Anhui","Hunan"),]) + 
  geom_line(aes(x=report_date,y=confirmed),col="red") + 
  geom_point(aes(x=report_date,y=confirmed),col="red") +
  geom_point(data=use_data_subset[use_data_subset$province %in% c("Hubei","Guangdong","Zheijiang","Henan","Anhui","Hunan"),],
             aes(x=report_date, y=confirmed),col="blue") +
  geom_line(data=use_data_subset[use_data_subset$province %in% c("Hubei","Guangdong","Zheijiang","Henan","Anhui","Hunan"),],
             aes(x=report_date, y=confirmed),col="blue") +
  scale_x_datetime(breaks="1 day") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  facet_wrap(~province,scales="free_y")



use_data_diff <- use_data_subset %>% 
  group_by(province) %>%
  mutate(prev=dplyr::lag(confirmed,n=1)) %>%
  mutate(prev=ifelse(raw_day==max(raw_day), NA, prev)) %>%
  mutate(diff=confirmed-prev) 

p_diff_per_day <- ggplot(use_data_diff) + 
    geom_bar(aes(x=raw_day,y=diff),stat="identity") + 
    facet_wrap(~province, scale="free_y",ncol=5) +
    theme_bw() +
    scale_x_date(limits=c(convert_date("19.01.2020"),convert_date(date_today)),
                 breaks="1 day") +
    xlab("Day of report") + ylab("New confirmed cases from previous day") +
    theme(axis.text.x=element_text(angle=45,hjust=1)) 
png("plots/reports_new_perday.png",width=12,height=12,res=300,units="in")
p_diff_per_day
dev.off()



