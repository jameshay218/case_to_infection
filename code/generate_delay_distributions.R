## Sliding confirmation date windows
## Generate a range of dates as a sliding window
date1 <- convert_date("01.01.2020")
date2 <- convert_date(date_today)
dates <- date1:date2

dates <- convert_date(dates)
range(kudos_dat$delay,na.rm=TRUE)

threshold <- 50

################################
## 1) PROBABILITIES GOING BACKWARD
## -- What's the probability that you developed symptoms on day t-x given confirmation on day t?

## Get number of new confirmations on each day
number_confirmed <- kudos_dat %>% filter(delay > 0) %>% group_by(reporting_date) %>% tally()

## Store generated counts and geometric probabilities for each day
all_dat_backward <- NULL

## For each day, go back in time until the number of new confirmations in that period
## is greater than threshold
first_date <- number_confirmed$reporting_date[which(cumsum(number_confirmed$n) > threshold)[1]]
use_dates <- dates[dates > first_date]
probs <- numeric(length(use_dates))
gamma_pars_backward <- NULL

model_probs_geometric_backward <- matrix(nrow=length(use_dates),ncol=101)
model_probs_gamma_backward <- matrix(nrow=length(use_dates),ncol=101)

## For each day of potential confirmation
for(i in seq_along(use_dates)) {
  counted <- 0
  start_date <- end_date <- use_dates[i]
  ## Go back in time and accumlate cases until >threshold cases
  while(counted < threshold & end_date > min(dates)){
   tmp_count <- number_confirmed %>% filter(reporting_date == end_date) %>% pull(n)
   counted <- counted + max(tmp_count,0)
   end_date <- end_date - 1
  }
  
  ## Get the data for this period
  tmp <- kudos_dat %>% filter(delay > 0 & reporting_date <= start_date &
                                   reporting_date > end_date ) %>% count(delay) %>% mutate(delay=as.numeric(delay))
  tmp <- tmp %>% mutate(start=use_dates[i])
  x <- kudos_dat %>% filter(delay > 0 & reporting_date <= start_date &
                              reporting_date > end_date) %>% pull(delay) %>% as.numeric
  ## VERY IMPORTANT
  ## We assume that it takes at least a day to get confirmed (ie. the ceiling of whatever
  ## randomly sampled confirmation delay time you got) 
  ## However, the fits are from a baseline of 0
  x_fit <- x - minimum_confirmation_delay
  
  ## Fit geometric distribution
  fit <- optim(0.15,fit_geometric, dat=x_fit,method="Brent",lower=0,upper=1)
  probs[i] <- fit$par
  
  ## Fit discretised gamma
  mean_start <- mean(x)
  var_start <- var(x)
  fit1 <- optim(par=c(5, 25), fn=fit_gamma_discrete,dat=x_fit)
  gamma_pars_backward[[i]] <- data.frame("gamma_mean_backward"=fit1$par[1],"gamma_var_backward"=fit1$par[2],
                                "date_confirmation"=use_dates[i], "n_used"=counted,
                                "direction"="backward")
  
  all_dat_backward <- bind_rows(tmp, all_dat_backward)
 
  model_probs_geometric_backward[i,] <- dgeom(0:100,probs[i])
  scale <- fit1$par[2]/fit1$par[1]
  shape <- fit1$par[1]/scale
  
  model_probs_gamma_backward[i,] <- ddgamma(0:100, scale=scale, shape=shape,log=FALSE)
}
all_dat_backward <- all_dat_backward %>% complete(delay, start, fill=list(n=0))

model_probs_geometric_backward <- reshape2::melt(model_probs_geometric_backward)
colnames(model_probs_geometric_backward) <- c("label", "delay","prob")
model_probs_geometric_backward$label <- use_dates[model_probs_geometric_backward$label]
model_probs_geometric_backward$label <- paste0("before ", model_probs_geometric_backward$label)
model_probs_geometric_backward$fit <- "Geometric"

model_probs_gamma_backward <- reshape2::melt(model_probs_gamma_backward)
colnames(model_probs_gamma_backward) <- c("label", "delay","prob")
model_probs_gamma_backward$delay <- model_probs_gamma_backward$delay
model_probs_gamma_backward$label <- use_dates[model_probs_gamma_backward$label]
model_probs_gamma_backward$label <- paste0("before ", model_probs_gamma_backward$label)
model_probs_gamma_backward$fit <- "Gamma"


all_dat_backward$label <- paste0("before ", all_dat_backward$start)

all_model_probs_backward <- bind_rows(model_probs_geometric_backward, model_probs_gamma_backward)

## For each day, go back in time day-by-day until at least 20 new confirmed cases have happened.
## Use the case confirmations in this window to generate a confirmation delay distribution
## for that window.

all_dat_backward <- all_dat_backward %>% group_by(label) %>% mutate(rel_n = n/sum(n))
p_sliding_delays_backward <- ggplot(all_dat_backward) + 
  geom_bar(aes(x=delay,y=rel_n),stat="identity") + 
  geom_line(data=all_model_probs_backward,aes(x=delay,y=prob, col=fit),size=1) +
  geom_vline(xintercept=1,linetype="dashed") +
  facet_wrap(~label) +
  coord_cartesian(xlim=c(0,40)) +
  theme_bw() +
  xlab("Delay from confirmation to symptom onset (days)") +
  ylab("Count") +
  ggtitle("Confirmation delay distribution from day of confirmation for each window (backward)") +
  theme(legend.position=c(0.9,0.1))
png("plots/delay_distribution_backward.png",height=6,width=8,units="in",res=300)
p_sliding_delays_backward
dev.off()

## What discretised gamma parameters should be used for each date of confirmation?
gamma_pars_dat_backward <- do.call("rbind", gamma_pars_backward)

################################
## 2) PROBABILITIES GOING FORWARD
## -- What's the probability that you have been confirmed by day t given symptom onset on day t-x?

## Get number of new onsets on each day
number_onsets <- kudos_dat %>% filter(delay > 0) %>% group_by(symptom_onset) %>% tally()

## Store generated counts and geometric probabilities for each day
all_dat_forward <- NULL

## For each day, go back in time until the number of new symptom onsets in that period
## is greater than threshold
first_date <- number_onsets$symptom_onset[which(cumsum(number_onsets$n) > threshold)[1]]
use_dates <- dates[dates > first_date]
probs <- numeric(length(use_dates))
gamma_pars_forward <- NULL

model_probs_geometric_forward <- matrix(nrow=length(use_dates),ncol=101)
model_probs_gamma_forward <- matrix(nrow=length(use_dates),ncol=101)

## For each day of potential confirmation
for(i in seq_along(use_dates)) {
  counted <- 0
  start_date <- end_date <- use_dates[i]
  ## Go back in time and accumlate cases until >threshold cases
  while(counted < threshold & end_date > min(dates)){
    tmp_count <- number_onsets %>% filter(symptom_onset == end_date) %>% pull(n)
    counted <- counted + max(tmp_count,0)
    end_date <- end_date - 1
  }
  
  ## Get the data for this period
  tmp <- kudos_dat %>% filter(delay > 0 & symptom_onset <= start_date &
                                symptom_onset > end_date ) %>% count(delay) %>% mutate(delay=as.numeric(delay))
  tmp <- tmp %>% mutate(start=use_dates[i])
  x <- kudos_dat %>% filter(delay > 0 & symptom_onset <= start_date &
                              symptom_onset > end_date) %>% pull(delay) %>% as.numeric
  ## VERY IMPORTANT
  ## We assume that it takes at least a day to get confirmed (ie. the ceiling of whatever
  ## randomly sampled confirmation delay time you got) 
  ## However, the fits are from a baseline of 0
  x_fit <- x - minimum_confirmation_delay
  
  ## Fit geometric distribution
  fit <- optim(0.15,fit_geometric, dat=x_fit,method="Brent",lower=0,upper=1)
  probs[i] <- fit$par
  
  ## Fit discretised gamma
  mean_start <- mean(x_fit)
  var_start <- var(x_fit)
  fit1 <- optim(par=c(5, 25), fn=fit_gamma_discrete,dat=x_fit)
  gamma_pars_forward[[i]] <- data.frame("gamma_mean_forward"=fit1$par[1],"gamma_var_forward"=fit1$par[2],
                                         "date_onset_symptoms"=use_dates[i], "n_used"=counted,
                                        "direction"="forward")
  
  all_dat_forward <- bind_rows(tmp, all_dat_forward)
  
  model_probs_geometric_forward[i,] <- dgeom(0:100,probs[i])
  scale <- fit1$par[2]/fit1$par[1]
  shape <- fit1$par[1]/scale
  
  model_probs_gamma_forward[i,] <- ddgamma(0:100, scale=scale, shape=shape,log=FALSE)
}
all_dat_forward <- all_dat_forward %>% complete(delay, start, fill=list(n=0))

model_probs_geometric_forward <- reshape2::melt(model_probs_geometric_forward)
colnames(model_probs_geometric_forward) <- c("label", "delay","prob")
model_probs_geometric_forward$label <- use_dates[model_probs_geometric_forward$label]
model_probs_geometric_forward$label <- paste0("before ", model_probs_geometric_forward$label)
model_probs_geometric_forward$fit <- "Geometric"

model_probs_gamma_forward <- reshape2::melt(model_probs_gamma_forward)
colnames(model_probs_gamma_forward) <- c("label", "delay","prob")
model_probs_gamma_forward$delay <- model_probs_gamma_forward$delay
model_probs_gamma_forward$label <- use_dates[model_probs_gamma_forward$label]
model_probs_gamma_forward$label <- paste0("before ", model_probs_gamma_forward$label)
model_probs_gamma_forward$fit <- "Gamma"


all_dat_forward$label <- paste0("before ", all_dat_forward$start)

all_model_probs_forward <- bind_rows(model_probs_geometric_forward, model_probs_gamma_forward)

## For each day, go back in time day-by-day until at least 20 new confirmed cases have happened.
## Use the case confirmations in this window to generate a confirmation delay distribution
## for that window.

all_dat_forward <- all_dat_forward %>% group_by(label) %>% mutate(rel_n = n/sum(n))
p_sliding_delays_forward <- ggplot(all_dat_forward) + 
  geom_bar(aes(x=delay,y=rel_n),stat="identity") + 
  geom_line(data=all_model_probs_forward,aes(x=delay,y=prob, col=fit),size=1) +
  geom_vline(xintercept=1,linetype="dashed") +
  coord_cartesian(xlim=c(0,40)) +
  facet_wrap(~label) +
  theme_bw() +
  xlab("Delay from symptom onset to confirmation (days)") +
  ylab("Count") +
  ggtitle("Confirmation delay distribution from day of symptom onset for each window (forward)") +
  theme(legend.position=c(0.9,0.1))
png("plots/delay_distribution_forward.png",height=6,width=8,units="in",res=300)
p_sliding_delays_forward
dev.off()

## What discretised gamma parameters should be used for each date of confirmation?
gamma_pars_dat_forward <- do.call("rbind", gamma_pars_forward)

