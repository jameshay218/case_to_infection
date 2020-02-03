// Stan model to fit geometric distribution to data
// data and fixed parameter values
data {
  int<lower=1> N; // number of samples
  int delay[N]; // delays
}

// model parameters
parameters {
  // ranges are ranges of uniform priors
  real<lower = 0., upper = 1.> p;
}

// observation model
// Stan does not implement the geometric distribution
// but the geometric distribution is a special case of the negative
// binomial distribution with the below parameterisation
model {
  for (n in 1:N) {
    delay[n] ~ neg_binomial(1, p/(1 - p));
  }
}
