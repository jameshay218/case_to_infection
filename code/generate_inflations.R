## How many symptom onsets do we have per day from confirmed cases?
symptom_observed <- sim_data_all %>% group_by(repeat_no, date_onset_symptoms) %>% tally()
colnames(symptom_observed)[2] <- "date"
symptom_observed$var <- "date_onset_symptoms"
## Confirm delay is time from symptom onset until now, as we want to find
## the probability that a case would have been confirmed by now from its onset
symptom_observed <- symptom_observed %>% mutate(confirm_delay=as.numeric(date_today-date))

## Need to inflate these, as only represent some proportion of actual symptom onsets based on how 
## long ago this was
if (use_geometric_confirmation_delay) {
  symptom_observed <- symptom_observed %>% left_join(confirm_probs_geometric)
} else {
  symptom_observed <- symptom_observed %>% left_join(confirm_probs_gamma)
}
## Because only some % of onsets from confirm_delay days ago will have been confirmed by now,
## we need to fill in the additional cases
symptom_observed <- symptom_observed %>% mutate(n_inflated=(rnbinom(n(), n, cumu_prob_confirm)))

## Now we have true inflated symptom onsets. For each symptom onset, on top of the infection times
## for the existing observed cases
symptom_unobserved <- symptom_observed %>% 
  select(repeat_no, date, n_inflated, confirm_delay) %>%
  group_by(date, repeat_no) %>% 
  uncount(n_inflated)

## Give each of these inflated symptom onsets an infection onset time
symptom_unobserved <- symptom_unobserved %>% left_join(used_weibull_pars)
symptom_unobserved <- symptom_unobserved %>% group_by(repeat_no) %>% mutate(augmented_id=1:n()) %>% ungroup()
symptom_unobserved <- symptom_unobserved %>% mutate(symp_delay = floor(rweibull(n(), alpha, sigma)),
                                                  date_infection = date - symp_delay,
                                                  total_delay = symp_delay + confirm_delay,
                                                  augmented=1,
                                                  individual=paste0("augmented_", repeat_no, "_", augmented_id)) %>%
  select(repeat_no, date, confirm_delay, symp_delay, total_delay, date_infection, individual, alpha, sigma, augmented) %>%
  ungroup()
colnames(symptom_unobserved)[2] <- "date_onset_symptoms"

## Now combine inflated and actual symptom onsets
## This is all cases with symptom onsets before today's date, inflated and observed
symptom_all <- sim_data_all %>% mutate(individual = as.character(individual)) %>% 
  ungroup() %>%
  bind_rows(symptom_unobserved)

symptom_all <- symptom_all %>% mutate(date_confirmation=ifelse(is.na(date_confirmation), 
                                                               date_onset_symptoms+confirm_delay, date_confirmation))
symptom_all <- symptom_all %>% mutate(date_confirmation = convert_date(date_confirmation))

## Tally infections per day with known symptom onset times
infections_with_symptoms <- symptom_all %>% group_by(repeat_no, date_infection) %>% tally()
## **********************
## VERY IMPORTANT
## this is augmenting infections starting from the first day that we could observe symptom onsets for
## i.e. yesterday, as minimum 1 day confirmation delay here
infections_with_symptoms <- infections_with_symptoms %>% mutate(symp_delay=as.numeric(date_today-date_infection)-minimum_confirmation_delay)
## **********************

## Now combine with symptom onset probs to find proportion of infections on each day
## that have not experienced symptoms by now. Then, get number of additional infections
## on each day that will not have experienced symptoms by now
infections_with_symptoms <- infections_with_symptoms %>% left_join(symptom_probs)
infections_with_symptoms <- infections_with_symptoms %>% mutate(n_inflated = rnbinom(n(), n, cumu_prob_symp))

## Spell out so that each row is an inflated case with an infection time 
## How many new infections roughly?
range_new_infections <- infections_with_symptoms %>%
  group_by(repeat_no) %>% 
  summarise(wow=sum(n_inflated)) %>% 
  pull(wow) %>% quantile(c(0.025,0.5,0.975))
## This many inflated infections (ie. upper and lower bound on unobserved infections)
## IN ADDITION to unobserved symptom onsets
print(range_new_infections)

## Expand out so 1 row per inflated infection
## This is inflated infections for each repeat
infections_with_symptoms <- infections_with_symptoms %>% group_by(repeat_no) %>% 
  mutate(augmented_id=1:n()) %>% ungroup()
infections_unobserved <- infections_with_symptoms %>% 
  select(repeat_no, date_infection, n_inflated, symp_delay, augmented_id) %>%
  group_by(date_infection, repeat_no) %>% 
  uncount(n_inflated) %>%
  mutate(augmented=1,
         individual=paste0("augmented_", repeat_no, "_", augmented_id)) %>%
  select(-augmented_id)

## Now merge unobserved infections with infections from those with symptom onset times
infections_all <- symptom_all %>% 
  select(repeat_no, individual, date_infection, symp_delay, augmented) %>%
  bind_rows(infections_unobserved)

## Get tallies to get bounds on infection and symptom onset numbers
## Stratified by whether case was augmented or not and total
final_infections_tally <- infections_all %>% 
  group_by(repeat_no, date_infection, augmented) %>%
  tally() %>% 
  pivot_wider(names_from=augmented,values_from=n, values_fill=list(n=0), names_prefix="inflated") %>%
  mutate(total=inflated0 + inflated1) %>% ungroup()
final_infections_tally$var <- "date_infections"
colnames(final_infections_tally)[2] <- "date"

## Stratified by whether case was augmented or not
final_symptom_onsets_tally<- symptom_all %>%
  group_by(repeat_no, date_onset_symptoms, augmented) %>%
  tally()  %>% 
  pivot_wider(names_from=augmented,values_from=n, values_fill=list(n=0), names_prefix="inflated") %>%
  mutate(total=inflated0 + inflated1) %>% ungroup()
final_symptom_onsets_tally$var <- "date_onset_symptoms"
colnames(final_symptom_onsets_tally)[2] <- "date"

## Now get bounds on these
final_all <- bind_rows(final_infections_tally, final_symptom_onsets_tally)
final_all <- final_all %>% complete(repeat_no, var, date, fill=list(inflated0=0,inflated1=0, total=0))

