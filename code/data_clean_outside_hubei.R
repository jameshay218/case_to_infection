## Get NON hubei data
github_dat_path <- "~/GitHub/nCoV2019/ncov_outside_hubei.csv"
other_dat <- read.csv(github_dat_path, stringsAsFactors=FALSE)
use_colnames <- c("age","sex","date_confirmation","date_onset_symptoms","date_admission_hospital")

## NON HUBEI DATA
#########
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


## Have a look at these variables over time
other_dat_tmp <- reshape2::melt(other_dat, id.vars=c("age","sex"))
other_dat_tmp <- other_dat_tmp %>% drop_na()
p_other_data <- ggplot(other_dat_tmp) + 
  geom_histogram(aes(x=value, fill=variable),stat="count") +
  facet_wrap(~variable, ncol=1,scales="free_y") +
  theme_bw() +
  scale_x_date(limits=c(convert_date("01.12.2019"),convert_date("31.01.2020")),
               breaks="2 day") +
  theme(axis.text.x=element_text(angle=45,hjust=1))
p_other_data

## Confirmation delay distribution from symptom onset
other_dat$confirmation_delay <- as.integer(other_dat$date_confirmation - other_dat$date_onset_symptoms)
other_dat$hospitalisation_delay <- as.integer(other_dat$date_admission_hospital - other_dat$date_onset_symptoms)

## One case looks like an error
other_dat1 <- other_dat[which(other_dat$hospitalisation_delay >= 0),]


p_other_hosp <- ggplot(other_dat1) + 
  geom_histogram(aes(x=hospitalisation_delay),stat="count") +
  scale_y_continuous(breaks=seq(0,nrow(other_dat),by=5)) +
  scale_x_continuous(breaks=seq(0,max(other_dat1$hospitalisation_delay + 10),by=1)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45,hjust=1))
p_other_hosp

other_dat2 <- other_dat[which(other_dat$confirmation_delay >= 0),]
p_other_confirm <- ggplot(other_dat2) + 
  geom_histogram(aes(x=confirmation_delay),stat="count") +
  scale_y_continuous(breaks=seq(0,nrow(other_dat),by=5)) +
  scale_x_continuous(breaks=seq(0,max(other_dat2$confirmation_delay + 10),by=1)) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45,hjust=1))
p_other_confirm

