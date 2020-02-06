library(incidence)

## Look at incidence of confirmed cases
## Get top provinces
inc_confirmed <- incidence(combined_dat_final$date_confirmation, groups=combined_dat_final$province)
fit_confirmed <- fit(inc_confirmed)
fit_confirmed
plot(inc_confirmed, fit=fit_confirmed)

paste0("Confirmed case reports doubling time: ", signif(fit_confirmed$info$doubling,3), " (95% CI; ",
       signif(fit_confirmed$info$doubling.conf[[1]],3),"-",signif(fit_confirmed$info$doubling.conf[[2]],3),")")

## Look at incidence of infections from simulated data
final_median <- final_all %>% group_by(var, date) %>% summarise(median_events=median(total))
final_median <- final_median %>% uncount(median_events)
## Get median infection onsets per day up to the 95% reported point
infections_median <- final_median %>% filter(var == "date_infections" & date <= thresholds[1]) 

inc_augmented <- incidence(infections_median$date)
fit_augmented <- fit(inc_augmented)
paste0("Augmented infection incidence doubling time: ", signif(fit_augmented$info$doubling,3), " (95% CI; ",
       signif(fit_augmented$info$doubling.conf[[1]],3),"-",signif(fit_augmented$info$doubling.conf[[2]],3),")")
plot(inc_augmented, fit=fit_augmented)


## REPEAT BY GROUP

## Look at incidence of confirmed cases
use_provinces <- c("Hubei","Zhejiang","Guangdong","Henan","Hunan","Jiangxi")
use_dat <- combined_dat_final %>% filter(province %in% use_provinces)
inc_confirmed <- incidence(use_dat$date_confirmation, groups=use_dat$province)
fit_confirmed <- fit(inc_confirmed)
fit_confirmed
plot(inc_confirmed, fit=fit_confirmed)

paste0("Confirmed case reports doubling time: ", signif(fit_confirmed$info$doubling,3), " (95% CI; ",
       signif(fit_confirmed$info$doubling.conf[[1]],3),"-",signif(fit_confirmed$info$doubling.conf[[2]],3),")")

## Look at incidence of infections from simulated data
final_use <- final_all_province %>% filter(province %in% use_provinces)
final_median <- final_use %>% group_by(var, date, province) %>% summarise(median_events=median(total))
final_median <- final_median %>% uncount(median_events)
## Get median infection onsets per day up to the 95% reported point
infections_median <- final_median %>% filter(var == "date_infections" & date <= thresholds[1]) 

inc_augmented <- incidence(infections_median$date, groups=infections_median$province)
fit_augmented <- fit(inc_augmented)
fit_augmented
paste0("Augmented infection incidence doubling time: ", signif(fit_augmented$info$doubling,3), " (95% CI; ",
       signif(fit_augmented$info$doubling.conf[[1]],3),"-",signif(fit_augmented$info$doubling.conf[[2]],3),")")
plot(inc_augmented, fit=fit_augmented)
