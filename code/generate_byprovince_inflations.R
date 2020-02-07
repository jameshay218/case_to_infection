#######################
## SPATIAL PLOTS
#######################
individual_key <- combined_dat_final[,c("individual","province","date_confirmation")]
sim_data_all <- as_tibble(sim_data_all)

sim_data_all_province <- right_join(individual_key, sim_data_all,by=c("individual","date_confirmation"))
sim_data_all_province <- sim_data_all_province %>% filter(!is.na(date_confirmation))
sim_data_all_province <- sim_data_all_province %>% select(individual, province, date_confirmation, 
                                      repeat_no, date_infection, date_onset_symptoms,
                                      alpha, sigma, symp_delay, confirm_delay, total_delay, augmented)

#start_dates <- sim_data_all_province %>% group_by(province,repeat_no) %>%
#  mutate(first_day=min(date_infection)) %>% ungroup()
## Days of first infection for each province
#p_start_days <- ggplot(start_dates) +
#  geom_histogram(aes(x=first_day)) + facet_wrap(~province,scales="free_y") +
#  theme_pubr()

## How many symptom onsets do we have per day from confirmed cases?
symptom_observed_province <- sim_data_all_province %>% group_by(repeat_no, date_onset_symptoms, province) %>% tally() %>% ungroup()
colnames(symptom_observed_province)[2] <- "date"
symptom_observed_province$var <- "date_onset_symptoms"
symptom_observed_province <- symptom_observed_province %>% mutate(confirm_delay=as.numeric(date_today-date))

## Need to inflate these, as only represent some proportion of actual symptom onsets based on how 
## long ago this was
if (use_geometric_confirmation_delay) {
  symptom_observed_province <- symptom_observed_province %>% left_join(confirm_probs_geometric)
} else {
  symptom_observed_province <- symptom_observed_province %>% left_join(confirm_probs_gamma)
}
## Because only some % of onsets from confirm_delay days ago will have been confirmed by now,
## we need to fill in the additional cases
symptom_observed_province <- symptom_observed_province %>% ungroup() %>% 
  group_by(province) %>% partition(cluster) 
symptom_observed_province <- symptom_observed_province %>%
  mutate(n_inflated=(rnbinom(n(), n, cumu_prob_confirm))) %>% collect()

## Now we have true inflated symptom onsets. For each symptom onset, on top of the infection times
## for the existing observed cases
symptom_unobserved_province <- symptom_observed_province %>% 
  select(repeat_no, date, n_inflated, confirm_delay, province) %>%
  group_by(date, repeat_no,province) %>% 
  uncount(n_inflated) %>% ungroup()

## Give each of these inflated symptom onsets an infection onset time
symptom_unobserved_province <- symptom_unobserved_province %>% left_join(used_weibull_pars)
symptom_unobserved_province <- symptom_unobserved_province %>% group_by(repeat_no) %>% mutate(augmented_id=1:n()) %>% ungroup()
symptom_unobserved_province <- symptom_unobserved_province %>% mutate(symp_delay = floor(rweibull(n(), alpha, sigma)),
                                                  date_infection = date - symp_delay,
                                                  total_delay = symp_delay + confirm_delay,
                                                  augmented=1,
                                                  individual=paste0("augmented_", repeat_no, "_", augmented_id)) %>%
  select(repeat_no, date, confirm_delay, symp_delay, total_delay, date_infection, 
         individual, alpha, sigma, province, augmented) %>%
  ungroup()
colnames(symptom_unobserved_province)[2] <- "date_onset_symptoms"

## Now combine inflated and actual symptom onsets
## This is all cases with symptom onsets before today's date, inflated and observed
symptom_all_province <- sim_data_all_province %>% mutate(individual = as.character(individual)) %>% 
  ungroup() %>%
  bind_rows(symptom_unobserved_province)
symptom_all_province <- symptom_all_province %>% mutate(date_confirmation=ifelse(is.na(date_confirmation), 
                                                               date_onset_symptoms+confirm_delay, date_confirmation))
symptom_all_province <- symptom_all_province %>% mutate(date_confirmation = convert_date(date_confirmation))

## Tally infections per day with known symptom onset times
infections_with_symptoms_province <- symptom_all_province %>% group_by(repeat_no, date_infection, province) %>% tally()
## **********************
## VERY IMPORTANT
## this is augmenting infections starting from the first day that we could observe symptom onsets for
## i.e. yesterday, as minimum 1 day confirmation delay here
infections_with_symptoms_province <- infections_with_symptoms_province %>% 
  mutate(symp_delay=as.numeric(date_today-date_infection)-minimum_confirmation_delay)
## **********************

## Now combine with symptom onset probs to find proportion of infections on each day
## that have not experienced symptoms by now. Then, get number of additional infections
## on each day that will not have experienced symptoms by now
infections_with_symptoms_province <- infections_with_symptoms_province %>% left_join(symptom_probs)
## Assign to cluster then run rbinoms
infections_with_symptoms_province <- infections_with_symptoms_province %>% group_by(province) %>% partition(cluster)
infections_with_symptoms_province <- infections_with_symptoms_province %>% mutate(n_inflated = rnbinom(n(), n, cumu_prob_symp)) %>% collect()
infections_with_symptoms_province <- infections_with_symptoms_province %>% group_by(repeat_no) %>% 
  mutate(augmented_id=1:n()) %>% ungroup()

## Expand out so 1 row per inflated infection
## This is inflated infections for each repeat
infections_unobserved_province <- infections_with_symptoms_province %>% 
  select(repeat_no, date_infection, n_inflated, symp_delay, province, augmented_id) %>%
  group_by(date_infection, repeat_no, province) %>% 
  uncount(n_inflated) %>%
  mutate(augmented=1,
         individual=paste0("augmented_", repeat_no, "_", augmented_id)) %>%
  select(-augmented_id)
## Now merge back with infections direct from symptom onset times
infections_all_province <- symptom_all_province %>% 
  select(repeat_no, individual, date_infection, symp_delay, augmented, province) %>%
  bind_rows(infections_unobserved_province)

## Get tallies to get bounds on infection and symptom onset numbers
## Stratified by whether case was augmented or not and total
final_infections_tally_province <- infections_all_province %>% 
  group_by(repeat_no, date_infection, augmented, province) %>%
  tally() %>% 
  pivot_wider(names_from=augmented,values_from=n, values_fill=list(n=0), names_prefix="inflated") %>%
  mutate(total=inflated0 + inflated1) %>% ungroup()
final_infections_tally_province$var <- "date_infections"
colnames(final_infections_tally_province)[2] <- "date"

## Stratified by whether case was augmented or not
final_symptom_onsets_tally_province <- symptom_all_province %>%
  group_by(repeat_no, date_onset_symptoms, augmented, province) %>%
  tally()  %>% 
  pivot_wider(names_from=augmented,values_from=n, values_fill=list(n=0), names_prefix="inflated") %>%
  mutate(total=inflated0 + inflated1) %>% ungroup()
final_symptom_onsets_tally_province$var <- "date_onset_symptoms"
colnames(final_symptom_onsets_tally_province)[2] <- "date"

## Now get bounds on these
final_all_province <- bind_rows(final_infections_tally_province, final_symptom_onsets_tally_province)
final_all_province <- final_all_province %>% complete(repeat_no, var, date, province, fill=list(inflated0=0,inflated1=0, total=0))

