## Pull latest cumulative incidence data from arcgis
url <- "https://docs.google.com/spreadsheets/d/1yZv9w9zRKwrGTaR-YzmAqMefw4wMlaXocejdxZaTs6w/htmlview?usp=sharing&sle=true"
available_dates <- sheets_sheets(url)

## Go through and pull each sheet
all_data <- NULL
for(date in available_dates) {
  tmp_dat <- read_sheet(ss=url,sheet=date)
  ## Some of the datetimes were entered incorrectly, so correct here
  if(date == "Jan25_12am") {
    tmp_dat <- tmp_dat %>% mutate(`Last Update`=as.POSIXct("2020-01-25 00:00:00 UTC",format="%Y-%m-%d %H:%M:%OS"))
  }
  if(date == "Jan25_10pm") {
    tmp_dat <- tmp_dat %>% mutate(`Last Update`=as.POSIXct("2020-01-25 22:00:00 UTC",format="%Y-%m-%d %H:%M:%OS"))
  }
  all_data <- bind_rows(all_data, tmp_dat)
}
use_data <- all_data

## More manageable variable names
colnames(use_data) <- c("province","country_region","updated",
                        "confirmed","deaths","recovered",
                        "suspected","demised","country",
                        "last_date")

## can't use data with no report date
use_data <- use_data %>% filter(!is.na(updated))
## Replace province with smallest admin unit with entry
use_data <- use_data %>% mutate(province=ifelse(is.na(province), country_region, province)) %>%
  mutate(province=ifelse(is.na(province), country, province))

## The other factors are in days, so can only use daily incidence.
## I reason that we should use the last report of each day that exists for
## each province
## Get the last report for each day (sorry this is such a rubbish solution)
subset_report_dates <- use_data %>% select(updated, province) %>%  
  unique %>% 
  arrange(province,updated) %>% 
  group_by(updated, province) %>%
  count() %>% ungroup() %>% 
  mutate(date=as.Date(updated, origin="1970-01-01")) %>%
  arrange(province, updated) %>%
  group_by(province,date) %>%
  mutate(number_reports=cumsum(n)) %>%
  filter(number_reports==max(number_reports)) %>% 
  arrange(province, updated) %>% 
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



p_cumu_confirmed_per_report <- ggplot(use_data) +
  geom_point(aes(x=as.Date(updated,origin="1970-01-01"),y=confirmed)) + 
  facet_wrap(~province,scales="free_y",ncol=5) +
  theme_bw() +
  scale_x_date(limits=c(convert_date("19.01.2020"),convert_date("31.01.2020")),
             breaks="1 day") +
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


