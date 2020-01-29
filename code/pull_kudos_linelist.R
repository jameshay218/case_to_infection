## Pulls the Kudos line list from google sheets
url <- "https://docs.google.com/spreadsheets/d/1jS24DjSPVWa4iuxuD4OAXrE3QeI8c9BC1hSlqr-NMiU/edit#gid=1187587451"
kudos_dat <- googlesheets4::read_sheet(url,sheet="Line-list",skip=1, na="NA")
kudos_dat$symptom_onset <- convert_date(unlist(lapply(kudos_dat$symptom_onset, convert_date)))
kudos_dat$hosp_visit_date <- convert_date(unlist(lapply(kudos_dat$hosp_visit_date, convert_date)))
kudos_dat$exposure_start <- convert_date(unlist(lapply(kudos_dat$exposure_start, convert_date)))
colnames(kudos_dat)[2] <- "reporting_date"

kudos_dat$reporting_date <- as.Date(kudos_dat$reporting_date)
kudos_dat <- kudos_dat %>% mutate(delay = reporting_date - symptom_onset)
kudos_dat <- kudos_dat %>% mutate(delay = ifelse(delay < 1, 1, delay))
