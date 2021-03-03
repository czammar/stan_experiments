model = """
data {
  int<lower=0> N;   // number of data items
  int<lower=0> K;   // number of predictors
  matrix[N, K] x;   // predictor matrix
  vector[N] y;      // outcome vector
}
parameters {
  real<lower=0> alpha;           // positive intercept
  vector[K] beta;       // coefficients for predictors
  real<lower=0> sigma;  // error scale
}
model {
  for (i in 1:K){beta[i] ~ normal(0, 1);} // priors
  y ~ normal(x * beta + alpha, sigma);  // likelihood
}
"""
