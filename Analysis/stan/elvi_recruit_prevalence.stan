data {
  int<lower=1> n_sites;
  int<lower=1> n_plots;
  int<lower=1> n_sites;

  
  vector[N] y;
}

parameters {
  real mu;
  real<lower=0> sigma;
}

model {
  y ~ normal(mu, sigma);
}

