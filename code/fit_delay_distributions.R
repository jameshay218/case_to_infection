
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
  ggtitle("Distribution of delays between symptom\n onset and confirmation") +
  theme_pubr()
p_other_confirm_fit

####################################
## OLD LINE LIST HOSPITALISATION DELAY DISTRIBUTION
####################################
## Assume that there is at least a 1 day delay to reporting, so < 1 day is set to 1
use_delays <- china_dat %>% select(hospitalisation_delay) %>% drop_na() %>% pull(hospitalisation_delay)
## Fit a geometric distribution to hospitalisation delay distribution
fit2 <- optim(c(0.1), fit_geometric, dat=use_delays,method="Brent",lower=0,upper=1)
times <- seq(0,25,by=1)
fit_line2 <- dgeom(times,prob=fit2$par)
fit_line_dat2 <- data.frame(x=times,y=fit_line2)

p_other_hosp_fit<- ggplot(china_dat) + 
  geom_histogram(aes(x=hospitalisation_delay,y=..density..),binwidth=1) +
  geom_line(data=fit_line_dat2, aes(x=x,y=y), col="red") +
  scale_x_continuous(breaks=seq(0,25,by=1)) +
  ggtitle("Distribution of delays between symptom\n onset and hospitalisation (not great fit)") +
  theme_pubr()
p_other_hosp_fit
## Fit isn't great for first day

####################################
## NEW LINE LIST HOSPITALISATION DELAY DISTRIBUTION
####################################
kudos_dat <- kudos_dat %>% mutate(hosp_delay = hosp_visit_date - symptom_onset)
## Assume that there is at least a 1 day delay to reporting, so < 1 day is set to 1
use_delays <- kudos_dat %>% select(hosp_delay) %>% drop_na() %>% pull(hosp_delay)
## Fit a geometric distribution to hospitalisation delay distribution
fit3 <- optim(c(0.1), fit_geometric, dat=use_delays,method="Brent",lower=0,upper=1)
times <- seq(0,25,by=1)
fit_line3 <- dgeom(times,prob=fit3$par)
fit_line_dat3 <- data.frame(x=times,y=fit_line3)

p_other_hosp_fit<- ggplot(kudos_dat) + 
  geom_histogram(aes(x=hosp_delay,y=..density..),binwidth=1) +
  geom_line(data=fit_line_dat2, aes(x=x,y=y), col="red") +
  scale_x_continuous(breaks=seq(0,25,by=1)) +
  ggtitle("Distribution of delays between symptom\n onset and hospitalisation (not great fit)") +
  theme_pubr()
p_other_hosp_fit


####################################
## KUDOS LINE LIST CONFIRMATION DELAY
####################################
## Fit a geometric distribution to the confirmation delay distribution
use_delays_kudos <- kudos_dat %>% select(delay) %>% drop_na() %>% pull(delay)
fit_kudos <- optim(c(0.1), fit_geometric, dat=use_delays_kudos-1,method="Brent",lower=0,upper=1)
fit_kudos_line <- dgeom(seq(0,max(kudos_dat$delay,na.rm=TRUE),by=1),prob=fit_kudos$par)
fit_line_kudos_dat <- data.frame(x=seq(1,max(kudos_dat$delay,na.rm=TRUE)+1,by=1),y=fit_kudos_line)

p_confirm_delay_kudos <- kudos_dat %>% select(delay) %>% drop_na() %>%
  ggplot() + 
  geom_histogram(aes(x=delay,y=..density..),binwidth=1,col="black") +
  geom_line(data=fit_line_kudos_dat, aes(x=x,y=y), col="red",size=1) +
  scale_x_continuous(breaks=seq(0,max(kudos_dat$delay,na.rm=TRUE),by=5),labels=seq(0,max(kudos_dat$delay,na.rm=TRUE),by=5)) +
  scale_y_continuous(expand=c(0,0),limits=c(0,0.15)) +
  geom_vline(xintercept=1,linetype="dashed") +
  ylab("Probability density") + xlab("Days since symptom onset") +
  ggtitle("Distribution of delays between symptom onset and confirmation\n Kudos line list data (geometric fit)") +
  theme_pubr()
p_confirm_delay_kudos


####################################
## KUDOS LINE LIST CONFIRMATION DELAY WITH GAMMA
####################################
## Fit a geometric distribution to the confirmation delay distribution
use_delays_kudos <- kudos_dat %>% select(delay) %>% drop_na() %>% pull(delay)
fit_kudos_gamma <- optim(c(5,25), fit_gamma_discrete, dat=use_delays_kudos-1)

fit_kudos_line_gamma <- dgamma_discrete_mean(seq(0,max(kudos_dat$delay,na.rm=TRUE),by=1),
                                             mean=fit_kudos_gamma$par[1],var=fit_kudos_gamma$par[2],use_log=FALSE)
fit_line_kudos_dat_gamma <- data.frame(x=seq(1,max(kudos_dat$delay,na.rm=TRUE)+1,by=1),y=fit_kudos_line_gamma)

p_confirm_delay_kudos_gamma <- kudos_dat %>% select(delay) %>% drop_na() %>%
  ggplot() + 
  geom_histogram(aes(x=delay,y=..density..),binwidth=1,col="black") +
  geom_line(data=fit_line_kudos_dat_gamma, aes(x=x,y=y), col="red",size=1) +
  scale_x_continuous(breaks=seq(0,max(kudos_dat$delay,na.rm=TRUE),by=5),labels=seq(0,max(kudos_dat$delay,na.rm=TRUE),by=5)) +
  scale_y_continuous(expand=c(0,0),limits=c(0,0.15)) +
  geom_vline(xintercept=1,linetype="dashed") +
  ylab("Probability density") + xlab("Days since symptom onset") +
  ggtitle("Distribution of delays between symptom onset and confirmation\n Kudos line list data (gamma fit)") +
  theme_pubr()
p_confirm_delay_kudos_gamma

