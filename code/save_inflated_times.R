inflated_infection_times <- sim_data_sum %>%filter(var == "date_infection" & repeat_no == 1)
inflated_infection_times <- unlist(Map(rep, inflated_infection_times$date, inflated_infection_times$n_inflated)) %>%
  as.data.frame
write_csv(inflated_infection_times, "augmented_data/inflated_infection_times.csv")

inflated_symptom_times <- sim_data_sum %>%filter(var == "date_onset_symptoms" & repeat_no == 1)
inflated_symptom_times <- unlist(Map(rep, inflated_symptom_times$date, inflated_symptom_times$n_inflated)) %>%
  as.data.frame
write_csv(inflated_symptom_times, "augmented_data/inflated_symptom_times.csv")
