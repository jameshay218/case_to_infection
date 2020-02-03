## Sliding confirmation date windows
## Generate a range of dates as a sliding window
date1 <- convert_date("01.01.2020")
date2 <- convert_date(date_today)-1
dates <- date1:date2

dates <- convert_date(dates)
range(kudos_dat$delay,na.rm=TRUE)

## Get number of new confirmations on each day
number_confirmed <- kudos_dat %>% filter(delay > 0) %>% group_by(reporting_date) %>% tally()

## Store generated counts and geometric probabilities for each day
all_dat <- NULL

## For each day, go back in time until the number of new confirmations in that period
## is greater than x=20
threshold <- 20
first_date <- number_confirmed$reporting_date[which(cumsum(number_confirmed$n) > threshold)[1]]
use_dates <- dates[dates >= first_date]
probs <- numeric(length(use_dates))
gamma_pars <- NULL

model_probs_geometric <- matrix(nrow=length(use_dates),ncol=41)
model_probs_gamma <- matrix(nrow=length(use_dates),ncol=41)

## For each day
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
  
  ## Fit geometric distribution
  fit <- optim(0.15,fit_geometric, dat=x-1,method="Brent",lower=0,upper=1)
  probs[i] <- fit$par
  
  ## Fit discretised gamma
  mean_start <- mean(x)
  var_start <- var(x)
  fit1 <- optim(par=c(5, 25), fn=fit_gamma_discrete,dat=x-1)
  gamma_pars[[i]] <- fit1$par
  
  all_dat <- bind_rows(tmp, all_dat)
 
  model_probs_geometric[i,] <- dgeom(0:40,probs[i])
  model_probs_gamma[i,] <- dgamma_mean_discrete(0:40, fit1$par[1], fit1$par[2],use_log=FALSE)
}
all_dat <- all_dat %>% complete(delay, start, fill=list(n=0))

model_probs_geometric <- reshape2::melt(model_probs_geometric)
colnames(model_probs_geometric) <- c("label", "delay","prob")
model_probs_geometric$label <- use_dates[model_probs_geometric$label]
model_probs_geometric$label <- paste0("before ", model_probs_geometric$label)
model_probs_geometric$fit <- "Geometric"

model_probs_gamma <- reshape2::melt(model_probs_gamma)
colnames(model_probs_gamma) <- c("label", "delay","prob")
model_probs_gamma$delay <- model_probs_gamma$delay
model_probs_gamma$label <- use_dates[model_probs_gamma$label]
model_probs_gamma$label <- paste0("before ", model_probs_gamma$label)
model_probs_gamma$fit <- "Gamma"


all_dat$label <- paste0("before ", all_dat$start)

all_model_probs <- bind_rows(model_probs_geometric, model_probs_gamma)

## For each day, go back in time day-by-day until at least 20 new confirmed cases have happened.
## Use the case confirmations in this window to generate a confirmation delay distribution
## for that window.

all_dat <- all_dat %>% group_by(label) %>% mutate(rel_n = n/sum(n))
p_sliding_delays <- ggplot(all_dat) + 
  geom_bar(aes(x=delay,y=rel_n),stat="identity") + 
  geom_line(data=all_model_probs,aes(x=delay,y=prob, col=fit),size=1) +
  geom_vline(xintercept=1,linetype="dashed") +
  facet_wrap(~label) +
  theme_bw() +
  xlab("Delay from symptom onset to confirmation (days)") +
  ylab("Count") +
  ggtitle("Delay distribution for each window")
png("plots/delay_distribution.png",height=6,width=8,units="in",res=300)
p_sliding_delays
dev.off()

