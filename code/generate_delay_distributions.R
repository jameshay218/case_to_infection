## Sliding confirmation date windows
## Generate a range of dates as a sliding window
date1 <- convert_date("01.01.2020")
date2 <- convert_date("28.01.2020")
dates <- date1:date2

dates <- convert_date(dates)
range(combined_dat$confirmation_delay,na.rm=TRUE)


## Get number of new confirmations on each day
number_confirmed <- combined_dat %>% filter(confirmation_delay > 0) %>% group_by(date_confirmation) %>% tally()

## Store generated counts and geometric probabilities for each day
all_dat <- NULL

## For each day, go back in time until the number of new confirmations in that period
## is greater than x
threshold <- 20
first_date <- number_confirmed$date_confirmation[which(cumsum(number_confirmed$n) > threshold)[1]]
use_dates <- dates[dates >= first_date]
probs <- numeric(length(use_dates))
model_probs <- matrix(nrow=length(use_dates),ncol=26)

for(i in seq_along(use_dates)) {
  counted <- 0
  start_date <- end_date <- use_dates[i]
  while(counted < threshold & end_date > min(dates)){
   tmp_count <- number_confirmed %>% filter(date_confirmation == end_date) %>% pull(n)
   counted <- counted + max(tmp_count,0)
   end_date <- end_date - 1
  }
  tmp <- combined_dat %>% filter(confirmation_delay > 0 & date_confirmation <= start_date &
                                 date_confirmation > end_date ) %>% count(confirmation_delay)
  x <- combined_dat %>% filter(confirmation_delay > 0 & date_confirmation <= start_date &
                            date_confirmation > end_date) %>% pull(confirmation_delay)
  fit <- optim(0.15,fit_geometric, dat=x,method="Brent",lower=0,upper=1)
  probs[i] <- fit$par
  tmp <- tmp %>% mutate(start=use_dates[i])
  all_dat <- bind_rows(tmp, all_dat)
  
  model_probs[i,] <- dgeom(0:25,probs[i])
}
plot(probs)
all_dat <- all_dat %>% complete(confirmation_delay, start, fill=list(n=0))

model_probs <- reshape2::melt(model_probs)
colnames(model_probs) <- c("label", "delay","prob")
model_probs$label <- use_dates[model_probs$label]
model_probs$label <- paste0("before ", model_probs$label)

all_dat$label <- paste0("before ", all_dat$start)

## For each day, go back in time day-by-day until at least 20 new confirmed cases have happened.
## Use the case confirmations in this window to generate a confirmation delay distribution
## for that window.

all_dat <- all_dat %>% group_by(label) %>% mutate(rel_n = n/sum(n))
p_sliding_delays <- ggplot(all_dat) + 
  geom_bar(aes(x=confirmation_delay,y=rel_n),stat="identity") + 
  geom_line(data=model_probs,aes(x=delay,y=prob),size=1,color="red") +
  geom_vline(xintercept=1,linetype="dashed") +
  facet_wrap(~label) +
  theme_bw() +
  xlab("Delay from symptom onset to confirmation (days)") +
  ylab("Count") +
  ggtitle("Delay distribution for each window")
