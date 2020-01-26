## Get hubei data
github_dat_path <- "~/GitHub/nCoV2019/ncov_hubei.csv"
hubei_dat <- read.csv(github_dat_path, stringsAsFactors=FALSE)


## For this, only need some of the variables
colnames(hubei_dat)
use_colnames <- c("age","sex","date_confirmation","date_onset_symptoms","date_admission_hospital")
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

## Have a look at these variables over time
hubei_dat_tmp <- reshape2::melt(hubei_dat, id.vars=c("age","sex"))
hubei_dat_tmp <- hubei_dat_tmp %>% drop_na()
p_hubei_data <- ggplot(hubei_dat_tmp) + 
  geom_histogram(aes(x=value, fill=variable),stat="count") +
  facet_wrap(~variable, ncol=1,scales="free_y") +
  theme_bw() +
  scale_x_date(limits=c(convert_date("01.12.2019"),convert_date("31.01.2020")),
               breaks="2 day") +
  theme(axis.text.x=element_text(angle=45,hjust=1))
p_hubei_data

## Clear outlier, culprit has wrong year
hubei_dat[which(hubei_dat$hospitalisation_delay < 0), "date_onset_symptoms"] <- c("2019-12-30","2020-01-18")

## Confirmation delay distribution from symptom onset
hubei_dat$confirmation_delay <- hubei_dat$date_confirmation - hubei_dat$date_onset_symptoms
hubei_dat$hospitalisation_delay <- as.integer(hubei_dat$date_admission_hospital - hubei_dat$date_onset_symptoms)

## One case looks like an error
hubei_dat <- hubei_dat[which(hubei_dat$hospitalisation_delay > 0),]

p_hubei_dist <- ggplot(hubei_dat) + 
  geom_histogram(aes(x=hospitalisation_delay),stat="count") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45,hjust=1))
p_hubei_dist
## The >30 day one is a bit dodgey
