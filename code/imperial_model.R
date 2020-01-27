## Shorter detection window is surely likely. If symptoms are onset before 
## travel, then case is surely more likely to be detected at incoming airport

## Pool of cases in Pop A
## All are infected suddenly
## Cases now have a probability of being detected in Pop A or in Pop B

## Prob of emerging in Pop A is mean time to detection, multipled by 
## the probability that they left within that time frame
## Every day, I have a probability of leaving
## If I leave, I stop trying to leave, otherwise I try again the next day
## So 1 - P(didn't leave)
## 1 - (1-P(leave))^days

number_cases <- 50
daily_passengers <- 3301
n_wuhan <- 19000000
detection_window <- 10
detection_prob <- 1

daily_prob <- daily_passengers/n_wuhan
prob_left <- daily_prob*detection_window

total_cases <- number_cases/(alt_prob*detection_prob)
total_cases







sum(p_x_k(seq(0,9,by=1),1,daily_prob))
p_x_k <- function(k, r, p){
  choose(k+r-1,k) * p^r * (1-p)^k
}

alt_prob <- 1 - (1-daily_prob)^detection_window

#total_cases1 <- number_cases/prob_left
#total_cases1


daily_prob +
  daily_prob*(1-daily_prob) +
  daily_prob*(1-daily_prob)^2 +
  daily_prob*(1-daily_prob)^3 +
  daily_prob*(1-daily_prob)^4 +
  daily_prob*(1-daily_prob)^5 +
  daily_prob*(1-daily_prob)^6 +
  daily_prob*(1-daily_prob)^7 +  
  daily_prob*(1-daily_prob)^8 +  
  daily_prob*(1-daily_prob)^9