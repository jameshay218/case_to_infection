## All the event delays (ie. incu + confirm delay)
## that would lead to an observation after a given delay
all_delays <- expand.grid(confirm=1:40, sympt=0:40)
all_delays$total_delay <- all_delays$confirm + all_delays$sympt
## Toy parameters for weibull dist on incubation period
alpha <- 2.379
sigma <- 6.613
days_past <- 41
## For each day from 41 days prior to current date, get the proportion of 
## infections observed given the possible associated delays, weighted by their
## probability of occurring
prop_seen <- numeric(days_past)
## For each day in the future
for(j in 1:days_past) {
  ## Get the event delay combinations that would lead
  ## to an infection being observed on this day
  tmp <- all_delays[all_delays$total_delay == j,]
  prob <- 0
  ## Probabilities for each of these combinations being observed on that day
  for(i in 1:nrow(tmp)){
    x <- tmp[i,]
    ## Weibull incubation period and geometric confirmation delay
    prob <- prob + dweibull(x$sympt, alpha, sigma) * dgeom(x$confirm-1, fit_kudos$par)
  }
  ## Proportion of past cases seen is therefore the sum of these probabilities
  prop_seen[j] <- prob
}
par(mar = c(5, 5, 5, 5))
total_days_considered <- length(prop_seen)
times <- times_unobserved <- seq(-total_days_considered, -1, 1)
#plot(times, rev(cumsum(prop_seen)), type = "l", 
#     xlab = "Date",
#     ylab = "Proportion of Infections Occurring That Day \nThat Have Been Reported By the Current Day")

