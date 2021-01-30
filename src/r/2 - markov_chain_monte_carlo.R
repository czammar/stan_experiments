############################################################
#
# Initial setup
#
############################################################

library(rstan)

c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

# To generate our samples we'll use R's pseudo random number
# generator which needs to be seeded to achieve reproducible
# results
set.seed(8675309)

# To ensure accurate results let's generate pretty large samples
N <- 100000

# To see how results scale with dimension we'll consider
# behavior one thorugh ten dimensions
Ds <- 1:10

idx <- rep(Ds, each=2)
plot_Ds <- sapply(1:length(idx), function(b) if(b %% 2 == 0) idx[b] + 0.5 else idx[b] - 0.5)

# Compute the MCMC summary statistics relevant to this exercise
fast_monitor <- function(sims, warmup) {
  dim_sims <- dim(sims)
  if (is.null(dim_sims)) {
    dim(sims) <- c(length(sims), 1, 1)
  } else if (length(dim_sims) == 2) {
    dim(sims) <- c(dim_sims, 1)
  } else if (length(dim_sims) > 3) {
    stop("'sims' has more than 3 dimensions")
  }
  
  parnames <- dimnames(sims)[[3]]
  if (is.null(parnames)) {
    parnames <- paste0("V", seq_len(dim(sims)[3]))
  }
  iter <- dim(sims)[1]
  chains <- dim(sims)[2]
  if (warmup > dim(sims)[1]) {
    stop("warmup is larger than the total number of iterations")
  }
  if (warmup >= 1) {
    sims <- sims[-seq_len(warmup), , , drop = FALSE]
  }
  
  out <- vector("list", length(parnames))
  out <- setNames(out, parnames)

  for (i in seq_along(out)) {
    sims_i <- sims[, , i]
    mean <- mean(sims_i)
    mcse_mean <- rstan:::mcse_mean(sims_i)
    ess <- round(rstan:::ess_rfun(sims_i))
    rhat <- rstan:::rhat_rfun(rstan:::split_chains(sims_i))
    out[[i]] <- c(mean, mcse_mean, ess, rhat)
  }
  
  out <- as.data.frame(do.call(rbind, out))
  colnames(out) <- c("mean", "se_mean", "n_eff", "Rhat")
  rownames(out) <- parnames
  
  out <- structure(out, chains = chains, iter = iter, warmup = warmup,
                   class = c("simsummary", "data.frame"))
  return(invisible(out))
}

############################################################
#
# How does the Random Walk Metropolis algorithm perform
# on a target distribution with a two-dimensional Gaussian 
# density function?
#
############################################################

# Target density
target_lpdf <- function(x) {
  - 0.5 * ( (x[1] - 1)**2 + (x[2] + 1)**2 ) - 0.5 * 2 * log(6.283185307179586)
}

# Tune proposal density
sigma <- 1.4

# A place to store our Markov chain
# D columns for the parameters and one extra column
# for the Metropolis acceptance probability
D <- 2 
mcmc_samples <- matrix(data = 0, nrow = N, ncol = D + 1)
  
# Randomly seed the initial state
mcmc_samples[1, 1:D] <- rnorm(D, 0, 3)
mcmc_samples[1, D + 1] <- 1
  
for (n in 2:N) {
  # Generate a proposal
  x0 <- mcmc_samples[n - 1, 1:D]
  xp <- rnorm(D, x0, sigma)
    
  # Compute acceptance probability
  accept_prob <- min(1, exp(target_lpdf(xp) - target_lpdf(x0)))
  mcmc_samples[n, D + 1] <- accept_prob
  
  # Apply Metropolis correction
  u = runif(1, 0, 1)
  if (accept_prob > u)
    mcmc_samples[n, 1:D] <- xp
  else
    mcmc_samples[n, 1:D] <- x0
}

# Compute MCMC estimator statistics, leaving
# out the first 100 samples as warmup
mcmc_stats <- fast_monitor(array(mcmc_samples, dim = c(N, 1, 3)), warmup=100)
rownames(mcmc_stats) <- c("x1", "x2", "accept_prob")  
mcmc_stats

# Plot convergence of MCMC estimators for each parameter
stride <- 200
M <- N / stride

iters <- stride * (1:(N / stride))

x1_mean <- rep(0, M)
x1_se <- rep(0, M)

x2_mean <- rep(0, M)
x2_se <- rep(0, M)

for (m in 1:M) {
  running_samples <- array(mcmc_samples[0:iters[m],], dim = c(iters[m], 1, 3))
  mcmc_stats <- fast_monitor(running_samples, warmup=100)
  
  x1_mean[m] <- mcmc_stats[1, "mean"]
  x1_se[m] <- mcmc_stats[1, "se_mean"]
  
  x2_mean[m] <- mcmc_stats[2, "mean"]
  x2_se[m] <- mcmc_stats[2, "se_mean"]
}

par(mfrow=c(1, 2))

plot(1, type="n", main="Mean of x1", 
     xlab="Iteration", xlim=c(0, N),
     ylab="Monte Carlo Estimator", ylim=c(0, 2))
  
polygon(c(iters, rev(iters)), c(x1_mean - 2 * x1_se, rev(x1_mean + 2 * x1_se)),
        col = c_light_highlight, border = NA)
lines(iters, x1_mean, col=c_dark, lwd=2)
abline(h=1, col="grey", lty="dashed", lw=2)


plot(1, type="n", main="Mean of x2", 
     xlab="Iteration", xlim=c(0, N),
     ylab="Monte Carlo Estimator", ylim=c(-2, 0))

polygon(c(iters, rev(iters)), c(x2_mean - 2 * x2_se, rev(x2_mean + 2 * x2_se)),
        col = c_light_highlight, border = NA)
lines(iters, x2_mean, col=c_dark, lwd=2)
abline(h=-1, col="grey", lty="dashed", lw=2)

############################################################
#
# How does the Random Walk Metropolis algorithm perform
# on a target distribution with a funnel density function?
#
############################################################

# Target density
target_lpdf <- function(x) {
  (- 0.5 * ( x[1]**2 + x[2]**2 + ((x[3] - x[1]) / exp(x[2]))**2 )
   - 0.5 * 3 * log(6.283185307179586) - 0.5 * x[2])
}

# Tune proposal density
sigma <- 1.4

# A place to store our Markov chain
# D columns for the parameters and one extra column
# for the Metropolis acceptance probability
D <- 3
mcmc_samples <- matrix(data = 0, nrow = N, ncol = D + 1)

# Randomly seed the initial state
mcmc_samples[1, 1:D] <- rnorm(D, 0, 3)
mcmc_samples[1, D + 1] <- 1

for (n in 2:N) {
  # Generate a proposal
  x0 <- mcmc_samples[n - 1, 1:D]
  xp <- rnorm(D, x0, sigma)
  
  # Compute acceptance probability
  accept_prob <- min(1, exp(target_lpdf(xp) - target_lpdf(x0)))
  mcmc_samples[n, D + 1] <- accept_prob
  
  # Apply Metropolis correction
  u = runif(1, 0, 1)
  if (accept_prob > u)
    mcmc_samples[n, 1:D] <- xp
  else
    mcmc_samples[n, 1:D] <- x0
}

# Compute MCMC estimator statistics, leaving
# out the first 100 samples as warmup
mcmc_stats <- fast_monitor(array(mcmc_samples, dim = c(N, 1, 4)), warmup=100)
rownames(mcmc_stats) <- c("x1", "x2", "x3", "accept_prob")  
mcmc_stats

# Plot convergence of MCMC estimators for the population parameters
stride <- 200
M <- N / stride

iters <- stride * (1:(N / stride))

mu_mean <- rep(0, M)
mu_se <- rep(0, M)

log_tau_mean <- rep(0, M)
log_tau_se <- rep(0, M)

for (m in 1:M) {
  running_samples <- array(mcmc_samples[0:iters[m],], dim = c(iters[m], 1, 4))
  mcmc_stats <- fast_monitor(running_samples, warmup=100)
  
  mu_mean[m] <- mcmc_stats[1, "mean"]
  mu_se[m] <- mcmc_stats[1, "se_mean"]
  
  log_tau_mean[m] <- mcmc_stats[2, "mean"]
  log_tau_se[m] <- mcmc_stats[2, "se_mean"]
}

par(mfrow=c(1, 2))

plot(1, type="n", main="Mean of mu", 
     xlab="Iteration", xlim=c(0, N),
     ylab="Monte Carlo Estimator", ylim=c(-1, 1))

polygon(c(iters, rev(iters)), c(mu_mean - 2 * mu_se, rev(mu_mean + 2 * mu_se)),
        col = c_light_highlight, border = NA)
lines(iters, mu_mean, col=c_dark, lwd=2)
abline(h=0, col="grey", lty="dashed", lw=2)


plot(1, type="n", main="Mean of log tau", 
     xlab="Iteration", xlim=c(0, N),
     ylab="Monte Carlo Estimator", ylim=c(-1, 1))

polygon(c(iters, rev(iters)), c(log_tau_mean - 2 * log_tau_se, rev(log_tau_mean + 2 * log_tau_se)),
        col = c_light_highlight, border = NA)
lines(iters, log_tau_mean, col=c_dark, lwd=2)
abline(h=0, col="grey", lty="dashed", lw=2)

############################################################
#
# How does the effective sample size of a Random Walk
# Metropolis Markov chain vary with the dimension of
# the target distribution?
#
############################################################

target_lpdf <- function(x) {
  -0.5 * sum(x**2) - 0.5 * length(x) * log(6.283185307179586)
}

############################################################
# First let's use a constant Markov transition
############################################################

accept_prob_means <- 0 * Ds
accept_prob_ses <- 0 * Ds
ave_eff_sample_sizes <- 0 * Ds

# Tune proposal density
sigma <- 1.4

for (D in Ds) {
  # A place to store our Markov chain
  # D columns for the parameters and one extra column
  # for the Metropolis acceptance probability
  mcmc_samples <- matrix(data = 0, nrow = N, ncol = D + 1)
  
  # Seeding the initial state with an exact sample
  # from the target distribution ensures that we
  # start in the typical set and avoid having to
  # worry about warmup.
  mcmc_samples[1, 1:D] <- rnorm(D, 0, 1)
  mcmc_samples[1, D + 1] <- 1
  
  for (n in 2:N) {
    # Generate a proposal
    x0 <- mcmc_samples[n - 1, 1:D]
    xp <- rnorm(D, x0, sigma)
    
    # Compute acceptance probability
    accept_prob <- min(1, exp(target_lpdf(xp) - target_lpdf(x0)))
    mcmc_samples[n, D + 1] <- accept_prob
    
    # Apply Metropolis correction
    u = runif(1, 0, 1)
    if (accept_prob > u)
      mcmc_samples[n, 1:D] <- xp
    else
      mcmc_samples[n, 1:D] <- x0
  }
  
  # Estimate average acceptance probability
  
  # Compute MCMC estimator statistics
  mcmc_stats <- fast_monitor(array(mcmc_samples, dim = c(N, 1, D + 1)), warmup=0)

  accept_prob_means[D] <- mcmc_stats[D + 1, "mean"]
  accept_prob_ses[D] <- mcmc_stats[D + 1, "se_mean"]
  
  # Estimate effective sample size
  ave_eff_sample_sizes[D] <- mean(mcmc_stats[1:D, "n_eff"])
}

par(mfrow=c(1, 2))

# Plot average acceptance probability verses dimension
pad_means <- do.call(cbind, lapply(idx, function(n) accept_prob_means[n]))
pad_ses <- do.call(cbind, lapply(idx, function(n) accept_prob_ses[n]))

plot(1, type="n", main="",
xlim=c(head(Ds, 1) - 0.5, tail(Ds, 1) + 0.5), xlab="Dimension",
ylim=c(0, 1), ylab="Average Acceptance Probability")

polygon(c(plot_Ds, rev(plot_Ds)), c(pad_means + 2 * pad_ses, rev(pad_means - 2 * pad_ses)),
col = c_light, border = NA)
lines(plot_Ds, pad_means, col=c_dark, lwd=2)

# Plot effective sample size per iteration verses dimension
pad_eff <- do.call(cbind, lapply(idx, function(n) ave_eff_sample_sizes[n]))

plot(1, type="n", main="", xlim=c(0, 10), xlab="Dimension",
ylim=c(0, 0.3), ylab="Average Effective Sample Size Per Iteration")

lines(plot_Ds, pad_eff / N, col=c_dark, lwd=2)

############################################################
# Now let's use an (approximately) optimally tuned Markov
# transition for each dimension
############################################################

accept_prob_means <- 0 * Ds
accept_prob_ses <- 0 * Ds
ave_eff_sample_sizes <- 0 * Ds

# Approximately optimal proposal tuning
opt_sigmas <- c(2.5, 1.75, 1.5, 1.2, 1.15, 1.0, 0.95, 0.85, 0.8, 0.75)

for (D in Ds) {
  # A place to store our Markov chain
  # D columns for the parameters and one extra column
  # for the Metropolis acceptance probability
  mcmc_samples <- matrix(data = 0, nrow = N, ncol = D + 1)
  
  # Seeding the initial state with an exact sample
  # from the target distribution ensures that we
  # start in the typical set and avoid having to
  # worry about warmup.
  mcmc_samples[1, 1:D] <- rnorm(D, 0, 1)
  mcmc_samples[1, D + 1] <- 1
  
  for (n in 2:N) {
    # Generate a proposal
    x0 <- mcmc_samples[n - 1, 1:D]
    xp <- rnorm(D, x0, opt_sigmas[D])
    
    # Compute acceptance probability
    accept_prob <- min(1, exp(target_lpdf(xp) - target_lpdf(x0)))
    mcmc_samples[n, D + 1] <- accept_prob
    
    # Apply Metropolis correction
    u = runif(1, 0, 1)
    if (accept_prob > u)
      mcmc_samples[n, 1:D] <- xp
    else
      mcmc_samples[n, 1:D] <- x0
  }
  
  # Estimate average acceptance probability
  
  # Compute MCMC estimator statistics
  mcmc_stats <- fast_monitor(array(mcmc_samples, dim = c(N, 1, D + 1)), warmup=0)
  
  accept_prob_means[D] <- mcmc_stats[D + 1, "mean"]
  accept_prob_ses[D] <- mcmc_stats[D + 1, "se_mean"]
  
  # Estimate effective sample size
  ave_eff_sample_sizes[D] <- mean(mcmc_stats[1:D, "n_eff"])
}

par(mfrow=c(1, 2))

# Plot average acceptance probability verses dimension
pad_means <- do.call(cbind, lapply(idx, function(n) accept_prob_means[n]))
pad_ses <- do.call(cbind, lapply(idx, function(n) accept_prob_ses[n]))

plot(1, type="n", main="",
     xlim=c(head(Ds, 1) - 0.5, tail(Ds, 1) + 0.5), xlab="Dimension",
     ylim=c(0, 1), ylab="Average Acceptance Probability")

polygon(c(plot_Ds, rev(plot_Ds)), c(pad_means + 2 * pad_ses, rev(pad_means - 2 * pad_ses)),
        col = c_light, border = NA)
lines(plot_Ds, pad_means, col=c_dark, lwd=2)

# Plot effective sample size per iteration verses dimension
pad_eff <- do.call(cbind, lapply(idx, function(n) ave_eff_sample_sizes[n]))

plot(1, type="n", main="", xlim=c(0, 10), xlab="Dimension",
     ylim=c(0, 0.3), ylab="Average Effective Sample Size Per Iteration")

lines(plot_Ds, pad_eff / N, col=c_dark, lwd=2)
