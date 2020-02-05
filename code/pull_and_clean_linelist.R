## Note this is called as part of "code/analysis_standalone.R"
save_my_plots <- TRUE

## Pulls up to date line list data from google sheets
##########################################
## 1. nCoV2019_2020_line_list_open
## Source: Xu et al. Epidemiological Data from the nCoV-2019 Outbreak: Early Descriptions from Publicly Available Data
## http://virological.org/t/epidemiological-data-from-the-ncov-2019-outbreak-early-descriptions-from-publicly-available-data/337
##########################################
url <- "https://docs.google.com/spreadsheets/d/1itaohdPiAeniCXNlntNztZ_oRvjh0HsGuJXUJWET008/edit?usp=sharing"
other_dat <- googlesheets4::read_sheet(url,sheet="outside_Hubei")
hubei_dat <- googlesheets4::read_sheet(url,sheet="Hubei")

## Look for confounding variable entries
print("Possible confounding date entries below. If any of these don't look like dates, add them to \"analysis_functions\": ")

## Range of allowable dates
valid_date_start <- convert_date("01.10.2019")
valid_date_end <- convert_date("01.01.2021")

valid_dates <- valid_date_start:valid_date_end
valid_dates <- convert_date(valid_dates)

test_colnames <- c("date_onset_symptoms", "date_admission_hospital", "date_confirmation",
                   "date_death_or_discharge")

all_failed_conversions <- NULL
all_outside_range <- NULL
for(i in seq_along(test_colnames)){
  message("Converting column: ", test_colnames[i])
  ## Outside hubei data
  unique_dates <- other_dat %>% pull(test_colnames[i])
  test1 <- convert_date(unique_dates)
  
  ## Get date entries outside allowable range
  outside_range <- test1[!(test1 %in% valid_dates) & !is.na(test1)]
  all_outside_range <- c(all_outside_range, outside_range)
  
  ## Date entries that could not be converted
  failed_conversions <- unique(unique_dates[intersect(which(!is.na(unique_dates)), which(is.na(test1)))])
  all_failed_conversions <- c(all_failed_conversions, failed_conversions)
  
  ## In Hubei data
  unique_dates <- hubei_dat %>% pull(test_colnames[i])
  test1 <- convert_date(unique_dates)
  
  ## Get date entries outside allowable range
  outside_range <- test1[!(test1 %in% valid_dates) & !is.na(test1)]
  all_outside_range <- c(all_outside_range, outside_range)
  
  ## Date entries that could not be converted
  failed_conversions <- unique(unique_dates[intersect(which(!is.na(unique_dates)), which(is.na(test1)))])
  all_failed_conversions <- c(all_failed_conversions, failed_conversions)
}
message(cat("Entries that could not be converted: ", unique(all_failed_conversions),sep="\n"))
message(cat("Entries that were outside range ", as.character(valid_date_start), "to", as.character(valid_date_end),": ", unique(all_outside_range),sep="\n"))

message(cat("Already excluded entries: ", confounding_dates, sep="\n"))

dput(unique(all_failed_conversions))
dput(unique(all_outside_range))
##########################################

############################
## 1. CLEAN DATES
############################
## For this, only need some of the variables
colnames(hubei_dat)
hubei_dat <- hubei_dat %>% select(c(use_colnames,"outcome","date_death_or_discharge"))
hubei_dat$age <- as.character(hubei_dat$age)
#########
## Clean up dates
#unique(hubei_dat$date_onset_symptoms)
#unique(hubei_dat$date_admission_hospital)
#unique(hubei_dat$date_confirmation)

hubei_dat$date_onset_symptoms <- clean_dates(hubei_dat$date_onset_symptoms)
hubei_dat$date_admission_hospital <- clean_dates(hubei_dat$date_admission_hospital)
hubei_dat$date_confirmation <- clean_dates(hubei_dat$date_confirmation)
hubei_dat$date_death_or_discharge <- clean_dates(hubei_dat$date_death_or_discharge)

## Convert to dates and remove those cases without known symptom onsets
hubei_dat$date_onset_symptoms <- convert_date(hubei_dat$date_onset_symptoms)
hubei_dat$date_admission_hospital <- convert_date(hubei_dat$date_admission_hospital)
hubei_dat$date_confirmation <- convert_date(hubei_dat$date_confirmation)
hubei_dat$date_death_or_discharge <- convert_date(hubei_dat$date_death_or_discharge)
hubei_dat$hubei <- 1

## NON HUBEI DATA
#########
other_dat$age <- as.character(other_dat$age)
## Clean up dates
other_dat <- other_dat %>% select(c(use_colnames,"date_death_or_discharge"))
#unique(other_dat$date_onset_symptoms)
#unique(other_dat$date_admission_hospital)
#unique(other_dat$date_confirmation)

other_dat$date_onset_symptoms <- clean_dates(other_dat$date_onset_symptoms)
other_dat$date_admission_hospital <- clean_dates(other_dat$date_admission_hospital)
other_dat$date_confirmation <- clean_dates(other_dat$date_confirmation)
other_dat$date_death_or_discharge <- clean_dates(other_dat$date_death_or_discharge)

## Convert to dates and remove those cases without known symptom onsets
other_dat$date_onset_symptoms <- convert_date(other_dat$date_onset_symptoms)
other_dat$date_admission_hospital <- convert_date(other_dat$date_admission_hospital)
other_dat$date_confirmation <- convert_date(other_dat$date_confirmation)
other_dat$date_death_or_discharge <- convert_date(other_dat$date_death_or_discharge)
other_dat$hubei <- 0
other_dat$outcome <- NA

combined_dat <- rbind(other_dat, hubei_dat)
#combined_dat <- combined_dat %>% select(c(use_colnames, "hubei"))
combined_dat$hubei <- as.factor(combined_dat$hubei)


## Replace province with smallest admin unit with entry
combined_dat <- combined_dat %>% mutate(province=ifelse(is.na(province), as.character(country), as.character(province)))


key_colnames <- key_colnames[key_colnames != "hubei"]
key_colnames <- c(key_colnames ,"hubei")

## Reorder factors by total number of confirmations
all_confirmations <- combined_dat %>% 
  select(country, date_confirmation) %>% 
  drop_na() %>%
  group_by(country) %>%
  tally() %>%
  arrange(-n)
factor_order <- all_confirmations %>% pull(country)
combined_dat$country <- factor(combined_dat$country, levels=factor_order)

china_dat <- combined_dat[combined_dat$country == "China",]

combined_dat_melted <- reshape2::melt(combined_dat, id.vars=c(key_colnames,"outcome"))


print("Generating line list data plots")
## Get more info on how many cases we have that are useable
p_missing_dist <- combined_dat_melted %>% 
  group_by(country, hubei, variable) %>% 
  filter(country == "China") %>%
  mutate(is_na = ifelse(is.na(value), "Missing", "Present")) %>%
  ggplot() + 
  geom_histogram(aes(x=is_na, fill=hubei),stat="count") +
  facet_grid(hubei~variable, scales="free_y") +
  xlab("Data missing?") +
  ggtitle("Missing data entries from crowdsources linelists \nfor China in/outside of Hubei") +
  scale_fill_manual(values=c("grey40","orange")) +
  theme_bw()

png("plots/linelist_missing_distribution.png",width=8,height=4,res=300,units="in")
p_missing_dist
dev.off()

combined_counts <- combined_dat_melted %>% group_by(province) %>% tally() %>% arrange(-n)
combined_dat_melted$province <- factor(combined_dat_melted$province, levels=combined_counts$province)
combined_dat_melted <- combined_dat_melted %>% drop_na()
combined_dat$province <- factor(combined_dat$province, levels=combined_counts$province)

## Have a look at these variables over time
p_other_data <- ggplot(combined_dat_melted[combined_dat_melted$variable %in% c("date_onset_symptoms","date_infections"),]) + 
  geom_histogram(aes(x=value, fill=hubei),stat="count") +
  facet_wrap(~variable, ncol=1,scales="free_y") +
  theme_bw() +
  scale_fill_manual(values=c("grey40","orange")) +
  scale_x_date(limits=c(convert_date("01.12.2019"),convert_date("31.01.2020")),
               breaks="7 day") +
  xlab("Date") + ylab("Count") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) + 
  facet_wrap(~country, scales="free_y",ncol=3)

png("plots/linelist_counts_country.png",width=10,height=6,res=300,units="in")
p_other_data
dev.off()
## Look by province
p_by_province <- ggplot(combined_dat_melted) + 
  geom_histogram(aes(x=value, fill=hubei),stat="count") +
  facet_wrap(~variable, ncol=1,scales="free_y") +
  theme_bw() +
  scale_fill_manual(values=c("grey40","orange")) +
  scale_x_date(limits=c(convert_date("01.12.2019"),convert_date("31.01.2020")),
               breaks="7 day") +
  xlab("Date") + ylab("Count") +
  theme(axis.text.x=element_text(angle=45,hjust=1),legend.position="none") + 
  facet_wrap(~province, scales="free_y",ncol=5)

png("plots/linelist_counts_province.png",width=10,height=10,res=300,units="in")
p_by_province
dev.off()

##########################################
## Get confirmation and hospitalisation delays
##########################################
combined_dat$confirmation_delay <- as.integer(combined_dat$date_confirmation - combined_dat$date_onset_symptoms)
combined_dat$hospitalisation_delay <- as.integer(combined_dat$date_admission_hospital - combined_dat$date_onset_symptoms)

## Look for any outliers
combined_dat %>% filter(combined_dat$confirmation_delay < 0)
combined_dat %>% filter(combined_dat$hospitalisation_delay < 0)

## Currently one person in China with negative hosp delay. Remove
combined_dat <- combined_dat %>% mutate(hospitalisation_delay = ifelse(hospitalisation_delay < 0, NA, hospitalisation_delay))
combined_dat <- combined_dat %>% mutate(confirmation_delay = ifelse(confirmation_delay < 0, NA, confirmation_delay))

## Have a look at delays by country, and hubei/outside hubei
## Confirmation delay
p_confirmation_delay_distribution <- ggplot(combined_dat) +
  geom_histogram(aes(x=confirmation_delay,fill=hubei),stat="count",binwidth=1) +
  scale_fill_manual(values=c("grey40","orange")) +
  ggtitle("Confirmation delay distribution from symptom onset") +
  theme_bw() +
  xlab("Days delay") + ylab("Count") +
  facet_wrap(~country,ncol=4,scales="free_y")

png("plots/linelist_conf_delay_dist.png",width=10,height=7,res=300,units="in")
p_confirmation_delay_distribution
dev.off()



p_confirmation_delay_china <- ggplot(combined_dat[combined_dat$country=="China",]) +
  geom_histogram(aes(x=confirmation_delay,fill=hubei),stat="count",binwidth=1) +
  scale_fill_manual(values=c("grey40","orange")) +
  ggtitle("Confirmation delay distribution China") +
  theme_bw() +
  xlab("Days delay from symptom onset") + ylab("Count") +
  facet_wrap(~hubei)

png("plots/linelist_conf_delay_china.png",width=8,height=4,res=300,units="in")
p_confirmation_delay_china
dev.off()


## Hospitalisation delay
p_hosp_delay_distribution <- ggplot(combined_dat) +
  geom_histogram(aes(x=hospitalisation_delay,fill=hubei),stat="count",binwidth=1) +
  scale_fill_manual(values=c("grey40","orange")) +
  theme_bw() +
  ggtitle("Hospitalisation delay distribution") +
  xlab("Days delay") + ylab("Count") +
  facet_wrap(~country,ncol=4,scales="free_y")

png("plots/linelist_hosp_delay_dist.png",width=10,height=7,res=300,units="in")
p_hosp_delay_distribution
dev.off()



p_hosp_delay_china <- ggplot(combined_dat[combined_dat$country=="China",]) +
  geom_histogram(aes(x=hospitalisation_delay,fill=hubei),stat="count",binwidth=1) +
  scale_fill_manual(values=c("grey40","orange")) +
  theme_bw() +
  xlab("Days delay") + ylab("Count") +
  ggtitle("Hospitalisation delay distribution China") +
  facet_wrap(~hubei)

png("plots/linelist_hosp_delay_china.png",width=8,height=4,res=300,units="in")
p_hosp_delay_china
dev.off()

##########################################
## Get death delay distribution
##########################################
combined_dat$death_delay <- as.integer(combined_dat$date_death_or_discharge - combined_dat$date_onset_symptoms)

## Look for any outliers
combined_dat %>% filter(combined_dat$death_delay < 0)

## Currently one person in China with negative hosp delay. Remove
combined_dat <- combined_dat %>% mutate(death_delay = ifelse(death_delay < 0, NA, death_delay))


## Death delay
p_death_delay_distribution <- ggplot(combined_dat) +
  geom_histogram(aes(x=death_delay,fill=hubei),stat="count",binwidth=1) +
  scale_fill_manual(values=c("grey40","orange")) +
  theme_bw() +
  ggtitle("Death delay distribution") +
  xlab("Days delay") + ylab("Count")

png("plots/linelist_death_delay_dist.png",width=10,height=7,res=300,units="in")
p_death_delay_distribution
dev.off()
