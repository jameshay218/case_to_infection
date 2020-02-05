## How many symptom onsets do we have per day from confirmed cases?
sim_symp_sum <- sim_data_all %>% group_by(repeat_no, date_onset_symptoms) %>% tally()
colnames(sim_symp_sum)[2] <- "date"
sim_symp_sum$var <- "date_onset_symptoms"
sim_symp_sum <- sim_symp_sum %>% mutate(confirm_delay=as.numeric(date_today-date))

## Need to inflate these, as only represent some proportion of actual symptom onsets based on how 
## long ago this was
sim_symp_sum <- sim_symp_sum %>% left_join(confirm_probs)
sim_symp_sum <- sim_symp_sum %>% mutate(n_inflated=(rnbinom(n(), n, cumu_prob_confirm)))

## Now we have true inflated symptom onsets. For each symptom onset, on top of the infection times
## for the existing observed cases
sim_symp_inflated <- sim_symp_sum %>% 
  select(repeat_no, date, n_inflated, confirm_delay) %>%
  group_by(date, repeat_no) %>% 
  uncount(n_inflated)

## Give each of these inflated symptom onsets an infection onset time
sim_symp_inflated <- sim_symp_inflated %>% left_join(used_weibull_pars)
sim_symp_inflated <- sim_symp_inflated %>% mutate(symp_delay = floor(rweibull(n(), alpha, sigma)),
                                                  date_infection = date - symp_delay,
                                                  total_delay = symp_delay + confirm_delay,
                                                  individual="augmented") %>%
  select(repeat_no, date, confirm_delay, symp_delay, total_delay, date_infection, individual, alpha, sigma) %>%
  ungroup()
colnames(sim_symp_inflated)[2] <- "date_onset_symptoms"

## Now combine inflated and actual symptom onsets
## This is all cases with symptom onsets before today's date, inflated and observed
sim_symp_full <- sim_data_all %>% mutate(individual = as.character(individual)) %>% 
  ungroup() %>%
  bind_rows(sim_symp_inflated)
sim_symp_full <- sim_symp_full %>% mutate(augmented = ifelse(individual=="augmented",1,0))

## Tally infections per day with known symptom onset times
sim_infections_symptoms <- sim_symp_full %>% group_by(repeat_no, date_infection) %>% tally()
sim_infections_symptoms <- sim_infections_symptoms %>% mutate(symp_delay=as.numeric(date_today-date_infection)-1)

## Now combine with symptom onset probs to find proportion of infections on each day
## that have not experienced symptoms by now. Then, get number of additional infections
## on each day that will not have experienced symptoms by now
sim_infections_symptoms <- sim_infections_symptoms %>% left_join(symptom_probs)
sim_infections_symptoms <- sim_infections_symptoms %>% mutate(n_inflated = rnbinom(n(), n, cumu_prob_symp))

## Spell out so that each row is an inflated case with an infection time 
## How many new infections roughly?
range_new_infections <- sim_infections_symptoms %>%
  group_by(repeat_no) %>% 
  summarise(wow=sum(n_inflated)) %>% 
  pull(wow) %>% quantile(c(0.025,0.5,0.975))
## This many inflated infections (ie. upper and lower bound on unobserved infections)
## IN ADDITION to unobserved symptom onsets
print(range_new_infections)

## Expand out so 1 row per inflated infection
## This is inflated infections for each repeat
sim_infections_symptoms_full <- sim_infections_symptoms %>% 
  select(repeat_no, date_infection, n_inflated, symp_delay) %>%
  group_by(date_infection, repeat_no) %>% 
  uncount(n_inflated) %>%
  mutate(individual="augmented_infection", augmented=1)

## Now merge back with infections direct from symptom onset times
final_infections <- sim_symp_full %>% 
  select(repeat_no, individual, date_infection, symp_delay, augmented) %>%
  bind_rows(sim_infections_symptoms_full)

## Total number of infections now much more reasonable
final_infections %>% group_by(repeat_no) %>% tally() %>% 
  pull(n) %>% quantile(c(0,0.025,0.5,0.975,1))

## Final number of symptom onsets
final_symptom_onsets <- sim_symp_full


## Get tallies to get bounds on infection and symptom onset numbers
## Stratified by whether case was augmented or not and total
final_infections_tally <- final_infections %>% 
  group_by(repeat_no, date_infection, augmented) %>%
  tally() %>% 
  pivot_wider(names_from=augmented,values_from=n, values_fill=list(n=0), names_prefix="inflated") %>%
  mutate(total=inflated0 + inflated1) %>% ungroup()
final_infections_tally$var <- "date_infections"
colnames(final_infections_tally)[2] <- "date"

## Stratified by whether case was augmented or not
final_symptom_onsets_tally<- final_symptom_onsets %>%
  group_by(repeat_no, date_onset_symptoms, augmented) %>%
  tally()  %>% 
  pivot_wider(names_from=augmented,values_from=n, values_fill=list(n=0), names_prefix="inflated") %>%
  mutate(total=inflated0 + inflated1) %>% ungroup()
final_symptom_onsets_tally$var <- "date_onset_symptoms"
colnames(final_symptom_onsets_tally)[2] <- "date"

## Now get bounds on these
final_all <- bind_rows(final_infections_tally, final_symptom_onsets_tally)
final_all <- final_all %>% complete(repeat_no, var, date, fill=list(inflated0=0,inflated1=0, total=0))

