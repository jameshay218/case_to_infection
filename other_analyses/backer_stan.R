## THIS IS ENTIRELY FROM Backer et al. medRxiv https://doi.org/10.1101/2020.01.27.20018986 

#########################################################################################
# Incubation period estimation for 2019-nCoV
# based on 34 confirmed cases outside of Wuhan in the period 20-23 January 2020
# with known travel history and symptom onset date
# 27 January 2020
# Jantien Backer (jantien.backer@rivm.nl)
#########################################################################################

library(tidyverse)
library(rstan)
library(loo)

rstan_options(auto_write = TRUE)
options(mc.cores = 10)

data <- read_tsv(file = "~/backer_s1_data.tsv")

tSymptomOnset <- as.integer((data$`Symptom onset` %>% as.Date(format = "%m/%d/%Y")) - as.Date("2019-12-31"))
tStartExposure <- as.integer((data$`Start exposure` %>% as.Date(format = "%m/%d/%Y")) - as.Date("2019-12-31"))
tEndExposure <- as.integer((data$`End exposure` %>% as.Date(format = "%m/%d/%Y")) - as.Date("2019-12-31"))

# exposure should ends before symptom onset
tEndExposure[tSymptomOnset < tEndExposure] <- tSymptomOnset[tSymptomOnset < tEndExposure]
# cases without maximum incubation period get start of exposure 21 days prior to first symptom onset
tStartExposure[is.na(tStartExposure)] <- min(tSymptomOnset) - 21

# input.data for stan analysis
input.data <- list(
  N = nrow(data),
  tStartExposure = tStartExposure,
  tEndExposure = tEndExposure,
  tSymptomOnset = tSymptomOnset)

# compile model

model <- stan(data = input.data, 
              chains = 0, 
              iter = 0,
              model_code = "
              data{
              int <lower = 1> N;
              vector[N] tStartExposure;
              vector[N] tEndExposure;
              vector[N] tSymptomOnset;
              }
              
              parameters{
              real<lower = 0> alphaInc; 	// Shape parameter of weibull distributed incubation period
              real<lower = 0> sigmaInc; 	// Scale parameter of weibull distributed incubation period
              vector<lower = 0, upper = 1>[N] uE;	// Uniform value for sampling between start and end exposure
              }
              
              transformed parameters{
              vector[N] tE; 	// infection moment
              tE = tStartExposure + uE .* (tEndExposure - tStartExposure);
              }
              
              model{
              // Contribution to likelihood of incubation period
              target += weibull_lpdf(tSymptomOnset -  tE  | alphaInc, sigmaInc);
              }
              
              generated quantities {
              // likelihood for calculation of looIC
              vector[N] log_lik;
              for (i in 1:N) {
              log_lik[i] = weibull_lpdf(tSymptomOnset[i] -  tE[i]  | alphaInc, sigmaInc);
              }
              }
              "
)

stanfit <- stan(fit = model, data = input.data, 
                init = "random",
                iter = 10000, chains = 8)


# check results and convergence
print(stanfit)

# modelfit
LL = extract_log_lik(stanfit, parameter_name = 'log_lik')
loo(LL)

# results
alpha <- rstan::extract(stanfit)$alphaInc
sigma <- rstan::extract(stanfit)$sigmaInc

# posterior median and 95%CI of mean
quantile(sigma*gamma(1+1/alpha), probs = c(0.025,0.5,0.975))

# posterior median and 95%CI of sd
quantile(sqrt(sigma^2*(gamma(1+2/alpha)-(gamma(1+1/alpha))^2)), probs = c(0.025,0.5,0.975))

# posterior median and 95%CI of percentiles
percentiles <- sapply(c(0.025, 0.05, 0.5, 0.95, 0.975, 0.99), function(p) quantile(qweibull(p = p, shape = alpha, scale = sigma), probs = c(0.025, 0.5, 0.975)))
colnames(percentiles) <- c(0.025, 0.05, 0.5, 0.95, 0.975, 0.99)

## JH ADDITIONS
weibull_draws <- data.frame(alpha=alpha, sigma=sigma)
write.csv(weibull_draws,"data/backer_weibull_draws.csv")
