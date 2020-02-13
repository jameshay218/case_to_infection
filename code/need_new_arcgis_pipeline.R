##'
##' Function to fit monotonically increasing spline. Returns a function that 
##' can be used to predict the spline for each date. 
##' 
##' 
##' @param dates dates of observations
##' @param obs the observations
##' @param df degrees of freedon. 
##' 
##' @return a function that takse in some number of dates and gives predictions on those
##' 
fit_ispline <- function (dates, obs, df=round(length(obs)/3)) {
  require(nnls)
  require(splines2)
  
  #first get the basis
  h <- iSpline(as.numeric(dates), df=df, intercept=T)
  
  
  #fit the nnls model to the data
  mdl <- nnls(h, obs)
  coefs <- coef(mdl)
  
  
  rc <- function(dates) {
    if(length(dates)==0) {return(NULL)}
    hnew <- predict(h, as.numeric(dates))
    return(hnew%*%coefs)
  }
  
  return(rc)
}


##'
##' Function to extract approximate epidemic curves
##' from the cumulative case data.
##'
##' @param cum_data a data frame with cumulative case data in oit
##' @param first_date the first date to infer
##' @param second_date the last date to infer...shold be iwthin data range
##'
##'
##' @return a data frame with roughly estimated incidence in it
##'
est_daily_incidence <- function (cum_data,
                                 first_date,
                                 last_date,
                                 na_to_zeros=FALSE) {
  if (na_to_zeros) {
    analyze <-   cum_data %>% replace(is.na(.), 0)
  } else {
    analyze <-   cum_data %>% drop_na(n)
  }
  
  
  ##Get the implied daily incidence for each province
  ##by fitting a monitonically increasing spline and then
  ##taking the difference (a little less sensitive to
  ##perturbations in reporting than taking raw difference).
  ##Making sure only to infer over trhe suport
  tmp_dt_seq <- seq(first_date, last_date, "days")
  incidence_data<- analyze %>% nest(-province) %>%
    mutate(cs=map(data, ~fit_ispline(dates=.$date, obs=.$n))) %>%
    mutate(Incidence=map2(cs,data, ~data.frame(Date=tmp_dt_seq[tmp_dt_seq>=min(.y$date) & tmp_dt_seq<=max(.y$date)],
                                               Incidence= diff(c(0, pmax(0,.x(tmp_dt_seq[tmp_dt_seq>=min(.y$date) & tmp_dt_seq<=max(.y$date)]))))))) %>%
    unnest(Incidence) %>% select(-data) %>% select(-cs)
  
  return(incidence_data)
  
  #####OLD VERSION 
  # tmp_dt_seq <- seq(first_date, last_date, "days")
  # incidence_data<- analyze %>% nest(-Province_State) %>%
  #   mutate(cs=map(data, ~splinefun(x=.$Update, y=.$Confirmed,
  #                                  method="hyman"))) %>%
  #   mutate(Incidence=map2(cs,data, ~data.frame(Date=tmp_dt_seq[tmp_dt_seq>=min(.y$Update) & tmp_dt_seq<=max(.y$Update)],
  #                                              Incidence= diff(c(0, pmax(0,.x(tmp_dt_seq[tmp_dt_seq>=min(.y$Update) & tmp_dt_seq<=max(.y$Update)]))))))) %>%
  #   unnest(Incidence) %>% select(-data) %>% select(-cs)
  
  return(incidence_data)
  
}


setwd("~/Documents/COVID-19/time_series/")

confirmed <- read_csv("time_series_2019-ncov-Confirmed.csv")
colnames(confirmed)[1:4] <- c("province","country_region","lat","long")
confirmed <- confirmed %>% pivot_longer(cols=-c("province","country_region","lat","long"), names_to="date",values_to="n")
confirmed$date <- as.POSIXct(confirmed$date,format="%m/%d/%Y %H:%M")

confirmed %>% filter(country_region == "Mainland China") %>%
ggplot() + geom_line(aes(x=date,y=n),stat="identity") + facet_wrap(~province,scales="free_y")

confirmed %>% filter(province=="Hubei") %>% ggplot() + geom_line(aes(x=date,y=n))

t_max <- lubridate::parse_date_time("2/08/20 00:00", c("%m/%d/%Y %I:%M %p", "%m/%d/%Y %H:%M", "%m/%d/%y %I:%M %p","%m/%d/%y %H:%M", "%Y-%m-%d %H:%M:%S"))
t_min <- lubridate::parse_date_time("1/22/20 00:00", c("%m/%d/%Y %I:%M %p", "%m/%d/%Y %H:%M", "%m/%d/%y %I:%M %p","%m/%d/%y %H:%M", "%Y-%m-%d %H:%M:%S"))
wow <- read_JHUCSSE_cases(t_max, FALSE)

prov <- "Zheijiang"
obs <- wow %>% arrange(Province_State, Update) %>% filter(Province_State==prov) %>% pull(Confirmed)
dates <- wow %>% arrange(Province_State, Update) %>% filter(Province_State==prov) %>% pull(Update)
#first get the basis
h <- iSpline(as.numeric(dates), df=df, intercept=T)
#fit the nnls model to the data
mdl <- nnls(h, obs)
coefs <- coef(mdl)
rc <- function(dates) {
  if(length(dates)==0) {return(NULL)}
  hnew <- predict(h, as.numeric(dates))
  return(hnew%*%coefs)
}
use_dates <- seq(min(dates), max(dates),by="days")

y <- rc(use_dates)

par(mfrow=c(2,1))
plot(obs~dates)
lines(y~use_dates,col="red")
plot(diff(obs)~dates[2:length(dates)])
lines(diff(y)~use_dates[2:length(use_dates)],col="red")


