library(tidyverse)
library(ggplot2)
library(patchwork)
dates <- 1:50
cases <- rpois(length(dates),floor((1:50)^1.7))
plot(cases~dates)

## Forward model
linelist <- tibble(n=cases,date_onset=dates) %>% uncount(n)
linelist <- linelist %>% mutate(delay=rgeom(n(), 0.15),
                                date_confirmation=date_onset + delay)
linelist <- linelist %>% filter(date_confirmation <= 50)


p1 <- linelist %>% ggplot() + geom_histogram(aes(x=delay),binwidth=1) + 
  facet_wrap(~date_onset,scales="free_y") + ggtitle("Forward")
p2 <- linelist %>% ggplot() + geom_histogram(aes(x=delay),binwidth=1) + 
  facet_wrap(~date_confirmation,scales="free_y") + ggtitle("Backward")

p1/p2

## Backward model
linelist_backward <- tibble(n=cases,date_confirmation=dates) %>% uncount(n)
linelist_backward <- linelist %>% mutate(delay=rgeom(n(), 0.15),
                                         date_onset=date_confirmation - delay)
linelist_backward <- linelist_backward %>% filter(date_confirmation <= 50 & date_confirmation > 0)


linelist_backward %>% ggplot() + geom_histogram(aes(x=delay),binwidth=1) + 
  facet_wrap(~date_onset,scales="free_y") + ggtitle("Forward")
linelist_backward %>% ggplot() + geom_histogram(aes(x=delay),binwidth=1) + 
  facet_wrap(~date_confirmation,scales="free_y") + ggtitle("Backward")


write_csv(combined_dat_final, "tmp.csv")


dat_ncov <- combined_dat_final$date_confirmation
real_dat <- tibble(date_confirmation=dat_ncov) %>% 
  mutate(delay=rgeom(n(), 0.15),
         date_onset=date_confirmation-delay)
real_dat <- real_dat %>% filter(date_confirmation > "2020-01-01")

real_dat %>% ggplot() + geom_histogram(aes(x=delay),binwidth=1) + 
  facet_wrap(~date_onset,scales="free_y") + ggtitle("Forward")
real_dat %>% ggplot() + geom_histogram(aes(x=delay),binwidth=1) + 
  facet_wrap(~date_confirmation,scales="free_y") + ggtitle("Backward")
