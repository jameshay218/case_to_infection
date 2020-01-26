Turning nCoV case reports into infection incidence
================
James Hay <jhay@hsph.harvard.edu>

## Introduction

Available line list data for SARS-nCoV mostly show the date of report
and not the date of infection onset. These reported cases therefore do
not directly represent underlying transmission dynamics, but rather the
day at which sick individuals become known to the healthcare system
(plus some delay in reporting). To understand the dynamics of the
outbreak, it is more informative to use the incidence curve for
infections rather than case reports. Here, I try to generate a more
realistic infection incidence curves for the aggregated case report data
from Wuhan.

## Method

Individuals undergo the following progression from infection:

![Timeline from infection to case reported used in the
model.](timeline.png)

I assume that on the day an individual seeks healthcare, they become
included in the case report number for that day. This assumption could
easily be relaxed, but it’s what I’ve assumed to start with. The
expected time of actual infection is therefore given by:

\(\text{expected day of infection} = \text{day you entered the hospital} - \text{expected delay from symptoms to healthcare} - \text{expected incubation period}\)

The number of cases reported on a given day is therefore the sum of
infections that occured on all preceding days multiplied by the
probability that those cases entered the hospital (and were therefore
reported) on that day of their infection.

Another way of thinking about it is that the number of infections on a
given day \(t\) is the sum of all cases that happened on all future days
\(i\), multipled by the probability of a case that started on day \(t\)
entering healthcare \(i-t\) days after infection.

\(X(t) = \sum_{i=t}^{t+t_{max}} C(i)*p(i-t)\)

where \(X(t)\) is the number of infections on day \(t\); \(C(i)\) is the
number of cases reported on day \(i\); and \(p(i-t)\) is the probability
of entering the hospital \(i-t\) days since the start of infection

I assume that individuals have 0 probability of entering the hospital
before symptom onset (\(\text{delay_onset}\)), and uniform daily
probability of entering the hospital in the following
\(\text{t_max_hosp}\) days. Again, this assumption could be relaxed to
use the actual delay distribution for symptom onset to hospitalisation,
but I’ve used a uniform distribution as proof of concept.

## Data

Data for the 41 confirmed cases pre- 01/01/2020 have already been
deaggregated into symptom onset times<sup>3</sup>. However, the cases
from 16/01/2020 (when the big spike of case reports started) are
aggregated by day. This analysis therefore applies to cases from
16/01/2020. Data were taken directly from a medarxiv
pre-print<sup>4</sup>. The same data (I believe) have also been
comendably compiled into line list data
[here](https://github.com/beoutbreakprepared/nCoV2019). Data were
accessed around midday 24/01/2020, so are likely now out of date.

``` r
wuhan_data <- read.csv("complete_timeline.csv",row.names=NULL)
wuhan_data$date <- as.Date(as.character(wuhan_data$date))
data <- wuhan_data[wuhan_data$date >= "20-01-16",]
print(wuhan_data)
```

    ##          date count
    ## 1  0019-12-01     1
    ## 2  0019-12-02     0
    ## 3  0019-12-03     0
    ## 4  0019-12-04     0
    ## 5  0019-12-05     0
    ## 6  0019-12-06     0
    ## 7  0019-12-07     0
    ## 8  0019-12-08     0
    ## 9  0019-12-09     0
    ## 10 0019-12-10     0
    ## 11 0019-12-11     3
    ## 12 0019-12-12     0
    ## 13 0019-12-13     0
    ## 14 0019-12-14     0
    ## 15 0019-12-15     0
    ## 16 0019-12-16     2
    ## 17 0019-12-17     0
    ## 18 0019-12-18     1
    ## 19 0019-12-19     1
    ## 20 0019-12-20     1
    ## 21 0019-12-21     5
    ## 22 0019-12-22     4
    ## 23 0019-12-23     3
    ## 24 0019-12-24     8
    ## 25 0019-12-25     1
    ## 26 0019-12-26     3
    ## 27 0019-12-27     2
    ## 28 0019-12-28     2
    ## 29 0019-12-29     0
    ## 30 0019-12-30     0
    ## 31 0019-12-31     0
    ## 32 0020-01-01     3
    ## 33 0020-01-02     1
    ## 34 0020-01-03     0
    ## 35 0020-01-04     0
    ## 36 0020-01-05     0
    ## 37 0020-01-06     0
    ## 38 0020-01-07     0
    ## 39 0020-01-08     0
    ## 40 0020-01-09     0
    ## 41 0020-01-10     0
    ## 42 0020-01-11     0
    ## 43 0020-01-12     0
    ## 44 0020-01-13     0
    ## 45 0020-01-14     0
    ## 46 0020-01-15     0
    ## 47 0020-01-16     4
    ## 48 0020-01-17    17
    ## 49 0020-01-18    59
    ## 50 0020-01-19    77
    ## 51 0020-01-20    60
    ## 52 0020-01-21    12
    ## 53 0020-01-22     0

## Parameters

``` r
delay_onset <- 5
t_max_hosp <- 8
prob_report <- 1
date_start <- as.Date("20-01-16")-(delay_onset + t_max_hosp)
```

## Results

This is the snipet of R code we’d change to use a different delay
distribution ie. `prob_report`.

``` r
#' Generate probability of seeking healthcare t_since_infection days after infection onset.
#'
#' @param t_since_infection days since infection started
#' @param delay_onset number of delays between infection and first day that you might seek healthcare
#' @param t_max_hosp maximum delay from symptom onset to reporting
#' @param overall_p overall probability that you will seek healthcare
prob_report <- function(t_since_infection, delay_onset=5, t_max_hosp=8,overall_p=1){
  prob_report_before_symptoms <- rep(0, delay_onset)
  prob_report_after_symptoms <- rep(overall_p/t_max_hosp, t_max_hosp)
  probs <- c(prob_report_before_symptoms, prob_report_after_symptoms)
  probs[t_since_infection]
}

## Store imputed data
imputed_data <- data.frame(date=seq(date_start, max(data$date),by=1),count=0)

## For each day
for(t in 1:nrow(imputed_data)){
  cases <- 0    
  date <- imputed_data[t,"date"]
  ## Get contribution of all future case reports to
  ## today's expected number of infections
  for(i in t:(t+(delay_onset + t_max_hosp-1))){
    check_date <- imputed_data[i,"date"]
    if(!is.na(check_date)){
      count <- wuhan_data[wuhan_data$date == check_date, "count"]
      prob <- prob_report(i-t+1, delay_onset, t_max_hosp)
      cases <- cases + count*prob
    }
  }
  imputed_data[t,"count"] <- cases
}
print(imputed_data)
```

    ##          date  count
    ## 1  0020-01-03  0.000
    ## 2  0020-01-04  0.500
    ## 3  0020-01-05  2.625
    ## 4  0020-01-06 10.000
    ## 5  0020-01-07 19.625
    ## 6  0020-01-08 27.125
    ## 7  0020-01-09 28.625
    ## 8  0020-01-10 28.625
    ## 9  0020-01-11 28.625
    ## 10 0020-01-12 28.125
    ## 11 0020-01-13 26.000
    ## 12 0020-01-14 18.625
    ## 13 0020-01-15  9.000
    ## 14 0020-01-16  1.500
    ## 15 0020-01-17  0.000
    ## 16 0020-01-18  0.000
    ## 17 0020-01-19  0.000
    ## 18 0020-01-20  0.000
    ## 19 0020-01-21  0.000
    ## 20 0020-01-22  0.000

## Incidence curve

``` r
## Replace aggregates with these counts
wuhan_data_old <- wuhan_data
times_in_both <- as.Date(intersect(wuhan_data$date,imputed_data$date),origin="1970-01-01")
wuhan_data[wuhan_data$date %in% times_in_both,"count"] <-
  round(imputed_data[imputed_data$date %in% times_in_both,"count"])

wuhan_data_old$ver <- "Old"
wuhan_data$ver <- "Raw"
overall_dat <- rbind(wuhan_data_old, wuhan_data)
ggplot(overall_dat) + 
  ggtitle("Data pre 01/01/2020 are date of symptom onset, data after are date of infection") +
  geom_bar(aes(x=date,y=count,fill=ver),position="dodge",stat="identity") +
  xlab("Date of report") +
  ylab("Infection incidence") +
  theme_bw()
```

![](README_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

## Comments

  - Modifying case data in this way might be useful as sensitivity
    analyses for models using the case reports as incidence (which is
    problematic if done naively)
  - We could use the actual [delay
    distribution](http://virological.org/t/epidemiological-data-from-the-ncov-2019-outbreak-early-descriptions-from-publicly-available-data/337)
    from the line list data rather than a uniform distribution over 8
    days
  - The reporting delay distribution has definitely changed over time.
    Again, this could be a function oof \(t\) rather than a single
    uniform distribution.

## Acknowledgements

Thanks Amy Dighe for checking that the logic here makes sense.

## References

1.  Cauchemez S, Fraser C, Van Kerkhove MD, Donnelly CA, Riley S,
    Rambaut A, et al. Middle East respiratory syndrome coronavirus:
    quantification of the extent of the epidemic, surveillance biases,
    and transmissibility. LANCET Infect Dis. 2014;14: 50–56.
    <doi:10.1016/S1473-3099(13)70304-9>
2.  Donnelly CA, Ghani AC, Leung GM, Hedley AJ, Fraser C, Riley S, et
    al. Epidemiological determinants of spread of causal agent of severe
    acute respiratory syndrome in Hong Kong. Lancet. 2003;361:
    1761–1766. <doi:10.1016/S0140-6736(03)13410-1>
3.  Huang et al. Clinical features of patients infected with 2019 novel
    coronavirus in Wuhan, China. Lancet 2020
    <doi:https://doi.org/10.1016/S0140-6736(20)30183-5>
4.  Read et al. Novel coronavirus 2019-nCoV: early estimation of
    epidemiological parameters and epidemic predictions. medRxiv
    <doi:https://doi.org/10.1101/2020.01.23.20018549>
