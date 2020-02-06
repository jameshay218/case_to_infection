####################################
## OLD LINE LIST CONFIRMATION DELAY DISTRIBUTION
####################################
china_dat <- combined_dat[combined_dat$country == "China",]
## Any instant confirmations?
china_dat <- china_dat %>% filter(confirmation_delay > 0 & confirmation_delay < 200)
## Assume that there is at least a 1 day delay to reporting, so < 1 day is set to 1
china_dat <- china_dat %>% mutate(confirmation_delay = ifelse(confirmation_delay < 1, 1, confirmation_delay))
use_delays <- china_dat %>% select(confirmation_delay) %>% filter(confirmation_delay < 200) %>% drop_na() %>% pull(confirmation_delay)
## Fit a geometric distribution to the confirmation delay distribution
fit1 <- optim(c(0.1), fit_geometric, dat=use_delays-1,method="Brent",lower=0,upper=1)
fit_line1 <- dgeom(seq(0,25,by=1),prob=fit1$par)
fit_line_dat1 <- data.frame(x=seq(1,26,by=1),y=fit_line1)

p_other_confirm_fit<- ggplot(china_dat) + 
  geom_histogram(aes(x=confirmation_delay,y=..density..),binwidth=1,col="black") +
  geom_line(data=fit_line_dat1, aes(x=x,y=y), col="red",size=1) +
  scale_x_continuous(breaks=seq(0,30,by=5),labels=seq(0,30,by=5)) +
  scale_y_continuous(expand=c(0,0),limits=c(0,0.2)) +
  geom_vline(xintercept=1,linetype="dashed") +
  ylab("Probability density") + xlab("Days since symptom onset") +
  ggtitle("Distribution of delays between symptom onset and confirmation") +
  theme_pubr()

####################################
## OLD LINE LIST HOSPITALISATION DELAY DISTRIBUTION
####################################
## Assume that there is at least a 1 day delay to reporting, so < 1 day is set to 1
use_delays <- china_dat %>% select(hospitalisation_delay) %>% 
  drop_na() %>% pull(hospitalisation_delay)
## Fit a geometric distribution to hospitalisation delay distribution
fit2 <- optim(c(0.1), fit_geometric, dat=use_delays,method="Brent",lower=0,upper=1)
times <- seq(0,25,by=1)
fit_line2 <- dgeom(times,prob=fit2$par)
fit_line_dat2 <- data.frame(x=times,y=fit_line2)

p_other_hosp_fit<- ggplot(china_dat) + 
  geom_histogram(aes(x=hospitalisation_delay,y=..density..),binwidth=1) +
  geom_line(data=fit_line_dat2, aes(x=x,y=y), col="red") +
  scale_x_continuous(breaks=seq(0,25,by=1)) +
  ggtitle("Distribution of delays between symptom onset and hospitalisation") +
  theme_pubr()
## Fit isn't great for first day

####################################
## NEW LINE LIST HOSPITALISATION DELAY DISTRIBUTION
####################################
## Assume that there is at least a 1 day delay to reporting, so < 1 day is set to 1
use_delays <- kudos_dat_china %>% select(hosp_delay) %>% drop_na() %>% pull(hosp_delay)
## Fit a geometric distribution to hospitalisation delay distribution
fit3 <- optim(c(0.1), fit_geometric, dat=use_delays,method="Brent",lower=0,upper=1)
times <- seq(0,25,by=1)
fit_line3 <- dgeom(times,prob=fit3$par)
fit_line_dat3 <- data.frame(x=times,y=fit_line3)

p_other_hosp_fit<- ggplot(kudos_dat_china) + 
  geom_histogram(aes(x=hosp_delay,y=..density..),binwidth=1) +
  geom_line(data=fit_line_dat3, aes(x=x,y=y), col="red") +
  scale_x_continuous(breaks=seq(0,25,by=1)) +
  ggtitle("Distribution of delays between symptom onset and hospitalisation") +
  theme_pubr()

####################################
## KUDOS LINE LIST CONFIRMATION DELAY
####################################
## Just China
## Fit a geometric distribution to the confirmation delay distribution
use_delays_kudos_global <- kudos_dat %>% select(delay) %>% drop_na() %>% pull(delay)
fit_kudos_global <- optim(c(0.1), fit_geometric, dat=use_delays_kudos_global-1,method="Brent",lower=0,upper=1)
fit_kudos_line_global <- dgeom(seq(0,max(kudos_dat$delay,na.rm=TRUE),by=1),prob=fit_kudos_global$par)
fit_line_kudos_dat_global <- data.frame(x=seq(1,max(kudos_dat$delay,na.rm=TRUE)+1,by=1),y=fit_kudos_line_global)

p_confirm_delay_kudos_global <- kudos_dat %>% select(delay) %>% drop_na() %>%
  ggplot() + 
  geom_histogram(aes(x=delay,y=..density..),binwidth=1,col="black") +
  geom_line(data=fit_line_kudos_dat_global, aes(x=x,y=y), col="red",size=1) +
  scale_x_continuous(breaks=seq(0,max(kudos_dat_china$delay,na.rm=TRUE),by=5),labels=seq(0,max(kudos_dat_china$delay,na.rm=TRUE),by=5)) +
  scale_y_continuous(expand=c(0,0),limits=c(0,0.15)) +
  geom_vline(xintercept=1,linetype="dashed") +
  ylab("Probability density") + xlab("Days since symptom onset") +
  ggtitle("Distribution of delays between symptom onset and confirmation\n Kudos line list data (geometric fit)") +
  theme_pubr()


## Fit a geometric distribution to the confirmation delay distribution
use_delays_kudos <- kudos_dat_china %>% select(delay) %>% drop_na() %>% pull(delay)
fit_kudos <- optim(c(0.1), fit_geometric, dat=use_delays_kudos-1,method="Brent",lower=0,upper=1)
fit_kudos_line <- dgeom(seq(0,max(kudos_dat_china$delay,na.rm=TRUE),by=1),prob=fit_kudos$par)
fit_line_kudos_dat <- data.frame(x=seq(1,max(kudos_dat_china$delay,na.rm=TRUE)+1,by=1),y=fit_kudos_line)

p_confirm_delay_kudos <- kudos_dat_china %>% select(delay) %>% drop_na() %>%
  ggplot() + 
  geom_histogram(aes(x=delay,y=..density..),binwidth=1,col="black") +
  geom_line(data=fit_line_kudos_dat, aes(x=x,y=y), col="red",size=1) +
  scale_x_continuous(breaks=seq(0,max(kudos_dat_china$delay,na.rm=TRUE),by=5),labels=seq(0,max(kudos_dat_china$delay,na.rm=TRUE),by=5)) +
  scale_y_continuous(expand=c(0,0),limits=c(0,0.15)) +
  geom_vline(xintercept=1,linetype="dashed") +
  ylab("Probability density") + xlab("Days since symptom onset") +
  ggtitle("Distribution of delays between symptom onset and confirmation\n Kudos line list data (geometric fit)") +
  theme_pubr()


####################################
## KUDOS LINE LIST CONFIRMATION DELAY WITH GAMMA
####################################
## Fit a geometric distribution to the confirmation delay distribution
use_delays_kudos <- kudos_dat_china %>% select(delay) %>% 
  drop_na() %>% pull(delay)
fit_kudos_gamma <- optim(c(5,25), fit_gamma_discrete, dat=use_delays_kudos-1)
fit_kudos_line_gamma <- dgamma_discrete_mean(seq(0,max(kudos_dat_china$delay,na.rm=TRUE),by=1),
                                             mean=fit_kudos_gamma$par[1],var=fit_kudos_gamma$par[2],use_log=FALSE)
fit_line_kudos_dat_gamma <- data.frame(x=seq(1,max(kudos_dat_china$delay,na.rm=TRUE)+1,by=1),y=fit_kudos_line_gamma)

p_confirm_delay_kudos_gamma <- kudos_dat_china %>% select(delay) %>% drop_na() %>%
  ggplot() + 
  geom_histogram(aes(x=delay,y=..density..),binwidth=1,col="black") +
  geom_line(data=fit_line_kudos_dat_gamma, aes(x=x,y=y), col="red",size=1) +
  scale_x_continuous(breaks=seq(0,max(kudos_dat_china$delay,na.rm=TRUE),by=5),labels=seq(0,max(kudos_dat_china$delay,na.rm=TRUE),by=5)) +
  scale_y_continuous(expand=c(0,0),limits=c(0,0.15)) +
  geom_vline(xintercept=1,linetype="dashed") +
  ylab("Probability density") + xlab("Days since symptom onset") +
  ggtitle("Distribution of delays between symptom onset and confirmation\n Kudos line list data (gamma fit)") +
  theme_pubr()





####################################
## KUDOS LINE LIST CONFIRMATION DELAY
####################################
## Fit a geometric distribution to the confirmation delay distribution
if(refit_p_confirm_delay) {
  use_delays_kudos <- kudos_dat_china %>% select(delay) %>% drop_na() %>% pull(delay)
  if(bayesian_p_confirm_delay) { # bayesian fit (posterior)
    confirmation_delay_model <- stan_model("code/geometric.stan")
    fit_kudos <- fit_geometric_stan(use_delays_kudos - 1, confirmation_delay_model) %>%
      extract(pars = "p")
    names(fit_kudos) <- "par"
  } else {
    # frequentist fit (point estimate)
    fit_kudos <- optim(c(0.1), fit_geometric, dat=use_delays_kudos-1,method="Brent",lower=0,upper=1)
  }
} else {
  ## or read from file
  if(bayesian_p_confirm_delay) {
    fit_kudos <- read.csv("data/p_confirm_delay_draws.csv") %>%
      as.list
    names(fit_kudos) <- "par"
  } else {
    fit_kudos <- list(par = as.numeric(read.csv("data/p_confirm_delay.csv")))
  }
}

plot_times <- seq(0,max(kudos_dat_china$delay,na.rm=TRUE))
predict_delay <- function(p_confirm_delay) {
  dgeom(plot_times, prob = p_confirm_delay)
} 

p_confirm_delay_kudos <- kudos_dat_china %>% select(delay) %>% drop_na() %>%
  ggplot() + 
  geom_histogram(aes(x=delay,y=..density..),binwidth=1,col="black") +
  scale_x_continuous(breaks=seq(0,max(kudos_dat$delay,na.rm=TRUE),by=5),labels=seq(0,max(kudos_dat$delay,na.rm=TRUE),by=5)) +
  scale_y_continuous(expand=c(0,0),limits=c(0,0.15)) +
  geom_vline(xintercept=1,linetype="dashed") +
  ylab("Probability density") + xlab("Days since symptom onset") +
  ggtitle("Distribution of delays between symptom onset\n and confirmation, Kudos line list data") +
  theme_pubr()
if(bayesian_p_confirm_delay) {
  fit_kudos_line <- vapply(fit_kudos$par, predict_delay, numeric(length(plot_times))) %>%
    apply(1, quantile, probs = c(0.025, 0.975))
  fit_line_kudos_dat <- data.frame(x=plot_times + 1,ymin=fit_kudos_line[1,],
                                   ymax = fit_kudos_line[2,])
  
  p_confirm_delay_kudos <- p_confirm_delay_kudos +
    geom_ribbon(data = fit_line_kudos_dat, aes(x = x, ymin = ymin, ymax = ymax), fill = "red")
} else {
  fit_kudos_line <- dgeom(seq(0,max(kudos_dat_china$delay,na.rm=TRUE),by=1),prob=fit_kudos$par)
  fit_line_kudos_dat <- data.frame(x=seq(1,max(kudos_dat_china$delay,na.rm=TRUE)+1,by=1),y=fit_kudos_line)
  
  p_confirm_delay_kudos <- p_confirm_delay_kudos +
    geom_line(data=fit_line_kudos_dat, aes(x=x,y=y), col="red",size=1)
}

p_confirm_delay_kudos


####################################
## SYMPTOM ONSET DISTRIBUTION
## UPDATED -- Now using the distribution from Weibull derived by
## Backer et al.
####################################
## 1000 draws from their posterior
n_samps <- 1000
times <- seq(0,25,by=0.1)
weibull_dists <- matrix(0, nrow=1000, ncol=length(times))
for(i in seq_len(n_samps)){
  pars <- weibull_stan_draws[i,]
  alpha <- pars$alpha
  sigma <- pars$sigma
  weibull_dists[i,] <- dweibull(times, alpha, sigma)
}
colnames(weibull_dists) <- times
weibull_dists_bounds <- as.data.frame(t(apply(weibull_dists, 2, function(x) quantile(x, c(0.025,0.5,0.975)))))
colnames(weibull_dists_bounds) <- c("lower","median","upper")
weibull_dists_bounds$times <- times

p_incubation <- ggplot(weibull_dists_bounds) + 
  geom_ribbon(aes(x=times, ymax=upper,ymin=lower),fill="grey70",alpha=0.4,col="black") + 
  geom_line(aes(x=times,y=median),size=1) +
  ylab("Probability density") +
  xlab("Days since onset of infection") +
  ggtitle("Incubation period distribution (Weibull, time from infection to symptoms)") +
  scale_y_continuous(limits=c(0,0.3),expand=c(0,0),breaks=seq(0,0.3,by=0.05)) +
  scale_x_continuous(expand=c(0,0)) +
  theme_pubr()


####################################
## DEATH DELAY DISTRIBUTION
####################################
## Fit a geometric distribution to the confirmation delay distribution
use_death_delays <- combined_dat %>% select(death_delay) %>% 
  drop_na() %>% pull(death_delay)
fit_deaths_gamma <- optim(c(5,25), fit_gamma_discrete, dat=use_death_delays)
fit_deaths_gamma_line <- dgamma_discrete_mean(seq(0,max(combined_dat$death_delay,na.rm=TRUE),by=1),
                                             mean=fit_deaths_gamma$par[1],var=fit_deaths_gamma$par[2],use_log=FALSE)
fit_deaths_gamma_dat <- data.frame(x=seq(0,max(combined_dat$death_delay,na.rm=TRUE),by=1),y=fit_deaths_gamma_line)

p_death_delay <- combined_dat %>% select(death_delay) %>% drop_na() %>%
  ggplot() + 
  geom_histogram(aes(x=death_delay,y=..density..),binwidth=1,col="black") +
  geom_line(data=fit_deaths_gamma_dat, aes(x=x,y=y), col="red",size=1) +
  scale_x_continuous(breaks=seq(0,max(kudos_dat_china$delay,na.rm=TRUE),by=5),labels=seq(0,max(kudos_dat_china$delay,na.rm=TRUE),by=5)) +
  #scale_y_continuous(expand=c(0,0),limits=c(0,0.15)) +
  ylab("Probability density") + xlab("Days from symptom onset to death") +
  theme_pubr()
