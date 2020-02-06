#' Vector of confounding date entries
confounding_dates <- c("none", "10.01.2020 - 22.01.2020", "pre 18.01.2020", "early january","",
                       "11.26.2020", "18.01.2020 - 23.01.2020", "not sure",-35067, 18660, 19025, ":",
                       "end of December","31.01,2020")


#' Checks the given vector of dates as strings and converts to NA if not useable
clean_dates <- function(dates){
  dates_new <- dates
  dates_new[dates_new %in% confounding_dates] <- NA
  dates_new
  
}
#' Converte date string to date format used in analysis
convert_date <- function(date1, format="%d.%m.%Y", origin="1970-01-01"){
  if(!(all(is.numeric(date1[!is.na(date1)])))) {
    ## Clean provided vector of dates
    date1 <- sapply(date1, function(x){
      clean <- x
      ## Keep as NA if needed
      if(!is.na(x)){
        ## Check for format dd.mm.yyyy or dd/mm/yyyy etc
        if(nchar(x) > 10){
          message(cat("Provided date \"", x, "\" must be less than 11 characters long\n"))
          clean <- NA
        } else {
          ## See if we can convert date. If it's something like "none", then will return NA
          clean <- tryCatch(as.Date(x, origin=origin, tryFormat=c(format)),
                            error = function(e) {
                              message(cat("Could not convert date \"", x, "\"\n"))
                              NA})
        }
      }
      clean
    }, simplify=FALSE)
    ## Seems redundant, but converting back to a vector needs going via numeric vector
    date1 <- as.numeric(unlist(date1))
  }
  as.Date(date1, origin=origin)
}


convert_datestring_to_date <- function(datestring){
  dates <- sapply(datestring, function(x){
    tmp <- strsplit(x, split="_")
    day <- as.numeric(str_extract(tmp[[1]][1], "[0-9]+"))
    month <- str_extract(tmp[[1]][1], "[aA-zZ]+")
    ## Extract time
    time <- str_extract(tmp[[1]][2], "[0-9]+")
    ## Only want first digit
    if(length(strsplit(time, "")[[1]]) == 3){
      time <- as.numeric(strsplit(time, "")[[1]][1])
    } else if(length(strsplit(time, "")[[1]]) == 4) {
      time <- as.numeric(paste(strsplit(time, "")[[1]][1:2], collapse= ""))
    } else {
      time <- as.numeric(time)
    }
    ## AM or PM?
    am <- str_extract(tmp[[1]][2], "[aA-zZ]+")
    if((am == "PM" | am == "pm") & time != 12) {
      time <- time + 12
    }
    if((am == "am" | am == "AM") & time == 12) {
      time <- 0
    }
    final_time <- make_datetime(2020, pmatch(month, month.name), day, time, 0)
    if(hour(final_time) == 0) final_time <- final_time - minutes(1)
    final_time
  },
  simplify=FALSE)
  do.call("c", dates)
}


## The updated ARC GIS dates mess with my pipeline, so for now manually plugging in the dates
sheet_name_key <- c("Feb02_9PM"="2020-02-02 21:00:00 UTC", 
                    "Feb02_745pm"="2020-02-02 19:45:00 UTC", 
                    "Feb02_5am"="2020-02-02 05:00:00 UTC",  
                    "Feb01_11pm"="2020-02-01 23:00:00 UTC",  
                    "Feb01_6pm"="2020-02-01 18:00:00 UTC",  
                    "Feb01_10am"="2020-02-01 10:00:00 UTC",  
                    "Jan31_7pm"="2020-01-31 19:00:00 UTC",  
                    "Jan31_2pm"="2020-01-31 14:00:00 UTC",  
                    "Jan30_930pm"="2020-01-30 21:30:00 UTC",  
                    "Jan30_11am"="2020-01-30 11:00:00 UTC",  
                    "Jan29_9pm"="2020-01-29 21:00:00 UTC",  
                    "Jan29_230pm"="2020-01-29 14:30:00 UTC",  
                    "Jan29_130pm"="2020-01-29 13:30:00 UTC",  
                    "Jan28_11pm"="2020-01-28 23:00:00 UTC",  
                    "Jan28_6pm"="2020-01-28 18:00:00 UTC",  
                    "Jan28_1pm"="2020-01-28 13:00:00 UTC", 
                    "Jan27_830pm"="2020-01-27 20:30:00 UTC", 
                    "Jan27_7pm"="2020-01-27 19:00:00 UTC",  
                    "Jan27_9am"="2020-01-27 09:00:00 UTC",  
                    "Jan26_11pm"="2020-01-26 23:00:00 UTC",  
                    "Jan26_11am"="2020-01-26 11:00:00 UTC", 
                    "Jan25_10pm"="2020-01-25 22:00:00 UTC",  
                    "Jan25_12pm"="2020-01-25 12:00:00 UTC",  
                    "Jan25_12am"="2020-01-24 23:59:59 UTC", 
                    "Jan24_12pm"="2020-01-24 12:00:00 UTC",  
                    "Jan24_12am"="2020-01-23 23:59:59 UTC",  
                    "Jan23_12pm"="2020-01-23 12:00:00 UTC",  
                    "Jan22_12pm"="2020-01-22 12:00:00 UTC",  
                    "Jan22_12am"="2020-01-21 23:59:59 UTC")