## Get hubei data
hubei_dat <- read.csv(hubei_data_path, stringsAsFactors=FALSE)

## For this, only need some of the variables
colnames(hubei_dat)
hubei_dat <- hubei_dat[,use_colnames]

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

## Get NON hubei data
other_dat <- read.csv(other_data_path, stringsAsFactors=FALSE)

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
other_dat$hubei <- 0

other_dat <- rbind(other_dat, hubei_dat)
other_dat$hubei <- as.factor(other_dat$hubei)
other_dat <- other_dat[other_dat$country == "China",]

key_colnames <- c(key_colnames ,"hubei")

print(paste0("Number of given confirmation dates (ie. max we can augment): ", nrow(other_dat[!is.na(other_dat$date_confirmation) | !is.na(other_dat$date_onset_symptoms),])))

## Have a look at these variables over time
other_dat_tmp <- reshape2::melt(other_dat, id.vars=key_colnames)
other_dat_tmp <- other_dat_tmp %>% drop_na()
p_other_data <- ggplot(other_dat_tmp) + 
  geom_histogram(aes(x=value, fill=hubei),stat="count") +
  facet_wrap(~variable, ncol=1,scales="free_y") +
  theme_bw() +
  scale_fill_manual(values=c("grey40","orange")) +
  scale_x_date(limits=c(convert_date("01.12.2019"),convert_date("31.01.2020")),
               breaks="2 day") +
  theme(axis.text.x=element_text(angle=45,hjust=1))
p_other_data


## Confirmation delay distribution from symptom onset
other_dat$confirmation_delay <- as.integer(other_dat$date_confirmation - other_dat$date_onset_symptoms)
other_dat$hospitalisation_delay <- as.integer(other_dat$date_admission_hospital - other_dat$date_onset_symptoms)

##########################################
## IMPORTANT, HANDLING OF OUTLIERS
##########################################

## Clear outliers, one culprit has wrong year, the other I'm not sure so will remove
other_dat[which(other_dat$hospitalisation_delay < 0), "date_onset_symptoms"] <- c("2019-12-30","2020-01-18")
other_dat <- other_dat[which(other_dat$hospitalisation_delay >= 0 | is.na(other_dat$hospitalisation_delay)),]

## Confirmation delay distribution from symptom onset
other_dat$confirmation_delay <- as.integer(other_dat$date_confirmation - other_dat$date_onset_symptoms)
other_dat$hospitalisation_delay <- as.integer(other_dat$date_admission_hospital - other_dat$date_onset_symptoms)

## One case looks like an error. Remove for now
other_dat1 <- other_dat[which(other_dat$hospitalisation_delay < 30),]
##########################################

p_other_hosp <- ggplot(other_dat1) + 
  geom_histogram(aes(x=hospitalisation_delay,fill=hubei),stat="count") +
  scale_y_continuous(breaks=seq(0,nrow(other_dat),by=5)) +
  scale_x_continuous(breaks=seq(0,max(other_dat1$hospitalisation_delay + 10),by=1)) +
  theme_bw() +
  scale_fill_manual(values=c("grey40","orange")) +
  theme(axis.text.x=element_text(angle=45,hjust=1))
p_other_hosp

other_dat2 <- other_dat[which(other_dat$confirmation_delay >= 0),]
p_other_confirm <- ggplot(other_dat2) + 
  geom_histogram(aes(x=confirmation_delay,fill=hubei),stat="count") +
  scale_y_continuous(breaks=seq(0,nrow(other_dat),by=5)) +
  scale_x_continuous(breaks=seq(0,max(other_dat2$confirmation_delay + 10),by=1)) +
  theme_bw() +
  scale_fill_manual(values=c("grey40","orange")) +
  theme(axis.text.x=element_text(angle=45,hjust=1))
p_other_confirm

