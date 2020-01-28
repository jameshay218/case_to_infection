############################
## 1. CLEAN DATES
############################
## For this, only need some of the variables
colnames(hubei_dat)
hubei_dat <- hubei_dat %>% select(use_colnames)
hubei_dat$age <- as.character(hubei_dat$age)
#########
## Clean up dates
unique(hubei_dat$date_onset_symptoms)
unique(hubei_dat$date_admission_hospital)
unique(hubei_dat$date_confirmation)

hubei_dat$date_onset_symptoms <- clean_dates(hubei_dat$date_onset_symptoms)
hubei_dat$date_admission_hospital <- clean_dates(hubei_dat$date_admission_hospital)
hubei_dat$date_confirmation <- clean_dates(hubei_dat$date_confirmation)

## Convert to dates and remove those cases without known symptom onsets
hubei_dat$date_onset_symptoms <- convert_date(hubei_dat$date_onset_symptoms)
hubei_dat$date_admission_hospital <- convert_date(hubei_dat$date_admission_hospital)
hubei_dat$date_confirmation <- convert_date(hubei_dat$date_confirmation)
hubei_dat$hubei <- 1

## NON HUBEI DATA
#########
other_dat$age <- as.character(other_dat$age)
## Clean up dates
other_dat <- other_dat[,use_colnames]
unique(other_dat$date_onset_symptoms)
unique(other_dat$date_admission_hospital)
unique(other_dat$date_confirmation)

other_dat$date_onset_symptoms <- clean_dates(other_dat$date_onset_symptoms)
other_dat$date_admission_hospital <- clean_dates(other_dat$date_admission_hospital)
other_dat$date_confirmation <- clean_dates(other_dat$date_confirmation)

## Convert to dates and remove those cases without known symptom onsets
other_dat$date_onset_symptoms <- convert_date(other_dat$date_onset_symptoms)
other_dat$date_admission_hospital <- convert_date(other_dat$date_admission_hospital)
other_dat$date_confirmation <- convert_date(other_dat$date_confirmation)
other_dat$hubei <- 0

combined_dat <- rbind(other_dat, hubei_dat)
combined_dat$hubei <- as.factor(combined_dat$hubei)

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
print(paste0("Number of given confirmation dates from China (ie. max we can augment): ", nrow(china_dat[!is.na(china_dat$date_confirmation) | !is.na(china_dat$date_onset_symptoms),])))

combined_dat_melted <- reshape2::melt(combined_dat, id.vars=key_colnames)

## Get more info on how many cases we have that are useable
p_missing_dist <- combined_dat_melted %>% 
  group_by(country, hubei, variable) %>% 
  filter(country == "China") %>%
  mutate(is_na = ifelse(is.na(value), "Missing", "Present")) %>%
  ggplot() + 
  geom_histogram(aes(x=is_na, fill=hubei),stat="count") +
  facet_grid(hubei~variable, scales="free_y") +
  xlab("Data missing?") +
  ggtitle("Number of missing data entries for China in/outside of Hubei") +
  scale_fill_manual(values=c("grey40","orange")) +
  theme_bw()
p_missing_dist

## Have a look at these variables over time
combined_dat_melted <- combined_dat_melted %>% drop_na()
p_other_data <- ggplot(combined_dat_melted) + 
  geom_histogram(aes(x=value, fill=hubei),stat="count") +
  facet_wrap(~variable, ncol=1,scales="free_y") +
  theme_bw() +
  scale_fill_manual(values=c("grey40","orange")) +
  scale_x_date(limits=c(convert_date("01.12.2019"),convert_date("31.01.2020")),
               breaks="7 day") +
  xlab("Date") + ylab("Count") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) + 
  facet_wrap(~country, scales="free_y")
p_other_data

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
  facet_wrap(~country)
p_confirmation_delay_distribution

p_confirmation_delay_china <- ggplot(combined_dat[combined_dat$country=="China",]) +
  geom_histogram(aes(x=confirmation_delay,fill=hubei),stat="count",binwidth=1) +
  scale_fill_manual(values=c("grey40","orange")) +
  ggtitle("Hospitalisation delay distribution") +
  theme_bw() +
  xlab("Days delay from symptom onset") + ylab("Count") +
  facet_wrap(~hubei)
p_confirmation_delay_china

## Hospitalisation delay
p_hosp_delay_distribution <- ggplot(combined_dat) +
  geom_histogram(aes(x=hospitalisation_delay,fill=hubei),stat="count",binwidth=1) +
  scale_fill_manual(values=c("grey40","orange")) +
  theme_bw() +
  xlab("Days delay") + ylab("Count") +
  facet_wrap(~country)
p_hosp_delay_distribution

p_hosp_delay_china <- ggplot(combined_dat[combined_dat$country=="China",]) +
  geom_histogram(aes(x=hospitalisation_delay,fill=hubei),stat="count",binwidth=1) +
  scale_fill_manual(values=c("grey40","orange")) +
  theme_bw() +
  xlab("Days delay") + ylab("Count") +
  facet_wrap(~hubei)
p_hosp_delay_china

