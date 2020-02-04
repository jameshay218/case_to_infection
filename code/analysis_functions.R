#' Colors for plots
blue_color <- "#0072B2"
orange_color <- "#D55E00"


#' Function for optim to fit poisson distribution
fit_poisson <- function(pars, dat){
  mean <- pars[1]
  -sum(dpois(dat, mean, TRUE))
}

################################
## GEOMETRIC FUNCTIONS
#' Function for optim to fit a geometric distribution
#' 
#' @param prob probability of success
#' @param dat is a vector of event observation times
#' @return negative sum log likelihood from dgeom
fit_geometric <- function(prob, dat){
  -sum(dgeom(x=dat, prob, log=TRUE))
}

#' fit geometric distribution using Stan
fit_geometric_stan <- function(delay_data, model) {
  
  # define data and constants
  data_list <- list(  
    N= length(delay_data),    # number of sampling times
    delay = delay_data)
  
  # fit the model to data
  chains <- 3 # number of NUTS chains to run in parallel
  set.seed(1) # set random number generator seed for reproducibility
  # fit model to data using No U-Turn sampler
  if(missing(model)) {
    fit <- stan(get_model_filename("geometric.stan"), data=data_list,  chains=3,verbose = TRUE)
  } else {
    fit <- sampling(model, data=data_list,  chains=3,verbose = TRUE)
  }
  fit
}


#################################
## GAMME FUNCTIONS
#' Gamma prob density from mean and variance
dgamma_mean <- function(x, mean, var, use_log=FALSE){
  scale <- var/mean
  shape <- mean/scale
  probs <- dgamma(x,shape=shape,scale=scale,log=use_log)
  probs
}

#' Gamma random generator from mean and variance
rgamma_mean <- function(n, mean, var){
  scale <- var/mean
  shape <- mean/scale
  rgamma(n,shape=shape,scale=scale)
}
#' Function for optim to fit a gamma distribution
#' 
#' @param pars vector, 1: gamma mean; 2: gamma variance
#' @param dat is a vector of event observation times
#' @return negative sum log likelihood from dgamma_mean
fit_gamma <- function(pars, dat){
  mean <- pars[1]
  var <- pars[2]
  -sum(dgamma_mean(dat, mean, var, TRUE))
}


#' Discretised gamma fit
dgamma_discrete_mean <- function(dat, mean, var, use_log=TRUE){
  scale <- var/mean
  shape <- mean/scale
  ddgamma(dat, shape=shape,scale=scale,log=use_log)
}

fit_gamma_discrete <- function(pars,dat){
  mean <- pars[1]
  var <- pars[2]
  scale <- var/mean
  shape <- mean/scale
  -sum(ddgamma(dat, shape=shape, scale=scale, log=TRUE))
  
}
