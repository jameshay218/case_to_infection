cluster <- new_cluster(4)
cluster_library(cluster, "tidyverse")

#######################
## SPATIAL PLOTS
#######################
individual_key <- combined_dat_final[,c("individual","province","date_confirmation")]
sim_data_all <- as_tibble(sim_data_all)

merged_data <- right_join(individual_key, sim_data_all,by=c("individual"))
merged_data <- merged_data %>% filter(!is.na(date_confirmation))
merged_data <- merged_data %>% select(individual, province, date_confirmation, 
                                      repeat_no, date_infection, date_onset_symptoms,
                                      alpha, sigma, symp_delay, confirm_delay, total_delay)

#start_dates <- merged_data %>% group_by(var, province,repeat_no) %>% 
#  mutate(first_day=min(date)) %>% ungroup() %>% 
#  group_by(var, province) %>%
#  filter(var=="date_infection") 

#ggplot(start_dates) +
#  geom_histogram(aes(x=first_day)) + facet_wrap(~province,scales="free_y")

## Get table of probabilities that events have happened after certain delays
if(!exists("all_probs_forward")) {
  all_probs_forward <- generate_forward_probabilities_dist(repeats, used_weibull_pars, 
                                                           fit_kudos$par,tmax=ceiling(max(sim_data_all$total_delay)))
}
confirm_probs <- all_probs_forward[[1]] %>% select(repeat_no, confirm_delay, cumu_prob_confirm)
symptom_probs <- all_probs_forward[[2]] %>% select(repeat_no, symp_delay, cumu_prob_symp)

## How many symptom onsets do we have per day from confirmed cases?
sim_symp_sum_province <- merged_data %>% group_by(repeat_no, date_onset_symptoms, province) %>% tally() %>% ungroup()
colnames(sim_symp_sum_province)[2] <- "date"
sim_symp_sum_province$var <- "date_onset_symptoms"
sim_symp_sum_province <- sim_symp_sum_province %>% mutate(confirm_delay=as.numeric(date_today-date))

## Need to inflate these, as only represent some proportion of actual symptom onsets based on how 
## long ago this was
sim_symp_sum_province <- sim_symp_sum_province %>% left_join(confirm_probs)
sim_symp_sum_province <- sim_symp_sum_province %>% ungroup() %>% group_by(province) %>% partition(cluster) 
sim_symp_sum_province <- sim_symp_sum_province %>%
  mutate(n_inflated=(rnbinom(n(), n, cumu_prob_confirm))) %>% collect()

## Now we have true inflated symptom onsets. For each symptom onset, on top of the infection times
## for the existing observed cases
sim_symp_inflated_province <- sim_symp_sum_province %>% 
  select(repeat_no, date, n_inflated, confirm_delay, province) %>%
  group_by(date, repeat_no,province) %>% 
  uncount(n_inflated) %>% ungroup()

## Give each of these inflated symptom onsets an infection onset time
sim_symp_inflated_province <- sim_symp_inflated_province %>% left_join(used_weibull_pars)
sim_symp_inflated_province <- sim_symp_inflated_province %>% mutate(symp_delay = floor(rweibull(n(), alpha, sigma)),
                                                  date_infection = date - symp_delay,
                                                  total_delay = symp_delay + confirm_delay,
                                                  individual="augmented") %>%
  select(repeat_no, date, confirm_delay, symp_delay, total_delay, date_infection, individual, alpha, sigma, province) %>%
  ungroup()
colnames(sim_symp_inflated_province)[2] <- "date_onset_symptoms"

## Now combine inflated and actual symptom onsets
## This is all cases with symptom onsets before today's date, inflated and observed
sim_symp_full_province <- merged_data %>% mutate(individual = as.character(individual)) %>% 
  ungroup() %>%
  bind_rows(sim_symp_inflated_province)
sim_symp_full_province <- sim_symp_full_province %>% mutate(augmented = ifelse(individual=="augmented",1,0))

rm(merged_data)

## Tally infections per day with known symptom onset times
sim_infections_symptoms_province <- sim_symp_full_province %>% group_by(repeat_no, date_infection, province) %>% tally()
sim_infections_symptoms_province <- sim_infections_symptoms_province %>% mutate(symp_delay=as.numeric(date_today-date_infection))

## Now combine with symptom onset probs to find proportion of infections on each day
## that have not experienced symptoms by now. Then, get number of additional infections
## on each day that will not have experienced symptoms by now
sim_infections_symptoms_province <- sim_infections_symptoms_province %>% left_join(symptom_probs)
## Assign to cluster then run rbinoms
sim_infections_symptoms_province <- sim_infections_symptoms_province %>% group_by(province) %>% partition(cluster)
sim_infections_symptoms_province <- sim_infections_symptoms_province %>% mutate(n_inflated = rnbinom(n(), n, cumu_prob_symp)) %>% collect()

## Spell out so that each row is an inflated case with an infection time 
## How many new infections roughly?
range_new_infections_province <- sim_infections_symptoms_province %>%
  group_by(repeat_no, province) %>% 
  summarise(wow=sum(n_inflated)) %>% 
  pull(wow) %>% quantile(c(0.025,0.5,0.975))
## This many inflated infections (ie. upper and lower bound on unobserved infections)
## IN ADDITION to unobserved symptom onsets
print(range_new_infections)

## Expand out so 1 row per inflated infection
## This is inflated infections for each repeat
sim_infections_symptoms_full_province <- sim_infections_symptoms_province %>% 
  select(repeat_no, date_infection, n_inflated, symp_delay, province) %>%
  group_by(date_infection, repeat_no, province) %>% 
  uncount(n_inflated) %>%
  mutate(individual="augmented_infection", augmented=1)

## Now merge back with infections direct from symptom onset times
final_infections_province <- sim_symp_full_province %>% 
  select(repeat_no, individual, date_infection, symp_delay, augmented, province) %>%
  bind_rows(sim_infections_symptoms_full_province)

## Final number of symptom onsets
final_symptom_onsets_province <- sim_symp_full_province

## Get tallies to get bounds on infection and symptom onset numbers
## Stratified by whether case was augmented or not and total
final_infections_tally_province <- final_infections_province %>% 
  group_by(repeat_no, date_infection, augmented, province) %>%
  tally() %>% 
  pivot_wider(names_from=augmented,values_from=n, values_fill=list(n=0), names_prefix="inflated") %>%
  mutate(total=inflated0 + inflated1) %>% ungroup()
final_infections_tally_province$var <- "date_infections"
colnames(final_infections_tally_province)[2] <- "date"

## Stratified by whether case was augmented or not
final_symptom_onsets_tally_province <- final_symptom_onsets_province %>%
  group_by(repeat_no, date_onset_symptoms, augmented, province) %>%
  tally()  %>% 
  pivot_wider(names_from=augmented,values_from=n, values_fill=list(n=0), names_prefix="inflated") %>%
  mutate(total=inflated0 + inflated1) %>% ungroup()
final_symptom_onsets_tally_province$var <- "date_onset_symptoms"
colnames(final_symptom_onsets_tally_province)[2] <- "date"

## Now get bounds on these
final_all_province <- bind_rows(final_infections_tally_province, final_symptom_onsets_tally_province)
final_all_province <- final_all_province %>% complete(repeat_no, var, date, province, fill=list(inflated0=0,inflated1=0, total=0))

