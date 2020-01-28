## Run this instead of pull_data.R to avoid using google.
## Need to have the data files saved locally as below.
hubei_data_path <- "~/Documents/nCoV2019/ncov_hubei.csv"
other_data_path <- "~/Documents/nCoV2019/ncov_outside_hubei.csv"

hubei_data_path <- "~/Documents/nCoV2019/dataset_archive/ncov_hubei_20200125 1906.csv"
other_data_path <- "~/Documents/nCoV2019/dataset_archive/ncov_outside_hubei_20200125 1848.csv"

## Get hubei data
hubei_dat <- read.csv(hubei_data_path, stringsAsFactors=FALSE)
## Get NON hubei data
other_dat <- read.csv(other_data_path, stringsAsFactors=FALSE)
