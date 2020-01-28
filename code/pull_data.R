source("code/analysis_functions.R")
library(tidyverse)
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

test_colnames <- c("date_onset_symptoms", "date_admission_hospital", "date_confirmation")

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
message(cat("Entries that could not be converted: ", unique(all_failed_conversions),sep=" "))
message(cat("Entries that were outside range ", as.character(valid_date_start), "to", as.character(valid_date_end),": ", unique(all_outside_range),sep=" "))

message(cat("Already excluded entries: ", confounding_dates, sep=" "))

dput(unique(all_failed_conversions))
dput(unique(all_outside_range))
##########################################


## Also want to use the Kudos line lists

