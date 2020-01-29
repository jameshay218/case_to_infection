## Toy parameters for weibull dist on incubation period
alpha <- 2.8
sigma <- 6

days_future <- 41
days_observed <- 31

## All the event delays (ie. incu + confirm delay)
## that would lead to an observation after a given delay
all_delays <- expand.grid(confirm=1:days_future, sympt=0:days_future)
all_delays$total_delay <- all_delays$confirm + all_delays$sympt

## For each day 0 to days_future - 1 days in the future
prop_seen <- numeric(days_future)

## For each day in the future
for(j in 1:days_future) {
  ## Get the event delay combinations that would lead
  ## to an infection being observed on this day
    tmp <- all_delays[all_delays$total_delay == j,]
    prob <- 0
    ## Probabilities for each of these combinations being observed on that day
    for(i in 1:nrow(tmp)){
      x <- tmp[i,]
      ## Weibull incubation and geometric confirmation delay
      prob <- prob + dweibull(x$sympt, alpha, sigma) * dgeom(x$confirm-1, 0.15)
    }
    ## Proportion of past cases seen is therefore the sum of these probabilities
    prop_seen[j] <- prob
}

## Now go forward in time for days_observed days to see how many existing infections
## we expect to have been observed
all_props <- matrix(0, ncol=days_observed,nrow=days_observed)
for(i in 1:days_observed){
  lag <- i
  capture <- 1:min(ncol(all_props) - i, length(prop_seen))
  all_props[i, capture+lag-1] <- prop_seen[capture]
}
plot(rowSums(all_props), type="l", 
     xlab="Days since start of outbreak",
     main="Proportion of infections that happened \non this day that could have been observed",
    ylab="Propn")
