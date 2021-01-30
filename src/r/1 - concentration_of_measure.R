############################################################
#
# Initial setup
#
############################################################

c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

# To facilitate the computation of Monte Carlo estimators let's
# define a _Welford accumulator_ that computes empirical means
# and variances of a sample in a single pass
welford_summary <- function(x) {
  summary = c(0, 0)
  for (n in 1:length(x)) {
    delta <- x[n] - summary[1]
    summary[1] <- summary[1] + delta / (n + 1)
    summary[2] <- summary[2] + delta * (x[n] - summary[1])
  }
  summary[2] <- summary[2] / (length(x) - 1)
  return(summary)
}

# We can then use the Welford accumulator output to compute the
# Monte Carlo estimator of a function and an estimate of its
# Monte Carlo Standard Error
compute_mc_stats <- function(x) {
  summary <- welford_summary(x)
  return(c(summary[1], sqrt(summary[2] / length(x))))
}

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


# Quantile probabilities that we'll use to quantify distributions
quant_probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

############################################################
#
# What is the volume of central rectangular box that spans
# [-1, +1] in each dimension relative to the volume of a
# box spanning [-3, +3] in each dimension?
#
############################################################

prob_means <- 0 * Ds
prob_ses <- 0 * Ds

for (D in Ds) {
  # Is the sampled point in the central interval?
  is_central_samples <- rep(0, N)
  
  for (n in 1:N) {
    # We start by assuming that the point will be
    # in the central interval
    is_central <- 1
    
    # Sample a new point one dimension at a time
    for (d in 1:D) {
      x_d <- runif(1, -3, 3)
      
      # If the component of the point in the current
      # dimension is not contained within the central
      # interval then set the flag to false
      if (-1 < x_d && x_d < 1)
      is_central <- is_central & 1
      else
      is_central <- is_central & 0
    }
    
    is_central_samples[n] <- is_central
  }
  
  # Estimate the relative volume as a probability
  s <- compute_mc_stats(is_central_samples)
  prob_means[D] <- s[1]
  prob_ses[D] <- s[2]
}

# Plot probabilities verses dimension
pad_means <- do.call(cbind, lapply(idx, function(n) prob_means[n]))
pad_ses <- do.call(cbind, lapply(idx, function(n) prob_ses[n]))

plot(1, type="n", main="",
xlim=c(head(Ds, 1) - 0.5, tail(Ds, 1) + 0.5), xlab="Dimension",
ylim=c(0, 0.5), ylab="Probability")

polygon(c(plot_Ds, rev(plot_Ds)), c(pad_means + 2 * pad_ses, rev(pad_means - 2 * pad_ses)),
col = c_light, border = NA)
lines(plot_Ds, pad_means, col=c_dark, lwd=2)

############################################################
#
# How much volume is in the neighborhood immediately outside
# a sphere, between a radius of 2 and 2.5, relative to the
# volume that lies in a neighborhood immediately inside that
# sphere, between a radius of 1.5 and 2?
#
############################################################

prob_inner_means <- 0 * Ds
prob_inner_ses <- 0 * Ds

prob_outer_means <- 0 * Ds
prob_outer_ses <- 0 * Ds

R <- 2
delta <- 0.5

for (D in Ds) {
  # Does the sampled point fall in the inside neighborhood?
  is_inner_samples <- rep(0, N)
  # Does the sampled point fall in the outside neighborhood?
  is_outer_samples <- rep(0, N)
  
  for (n in 1:N) {
    # Sample a new point
    x <- runif(D, -3, 3)
    
    # Compute distance from origin
    r <- sqrt(sum(x**2))
    
    # Check if point falls in the inside neighborhood
    if (R - delta < r && r < R)
    is_inner_samples[n] <- 1;
    
    # Check if point falls in the outside neighborhood
    if (R < r && r < R + delta)
    is_outer_samples[n] <- 1;
  }
  
  # Estimate the relative volumes as probabilies
  s1 <- compute_mc_stats(is_inner_samples)
  prob_inner_means[D] <- s1[1]
  prob_inner_ses[D] <- s1[2]
  
  s2 <- compute_mc_stats(is_outer_samples)
  prob_outer_means[D] <- s2[1]
  prob_outer_ses[D] <- s2[2]
}

# Plot probabilities verses dimension
pad_inner_means <- do.call(cbind, lapply(idx, function(n) prob_inner_means[n]))
pad_inner_ses <- do.call(cbind, lapply(idx, function(n) prob_inner_ses[n]))

pad_outer_means <- do.call(cbind, lapply(idx, function(n) prob_outer_means[n]))
pad_outer_ses <- do.call(cbind, lapply(idx, function(n) prob_outer_ses[n]))

plot(1, type="n", main="",
xlim=c(head(Ds, 1) - 0.5, tail(Ds, 1) + 0.5), xlab="Dimension",
ylim=c(0, 0.25), ylab="Probability")

polygon(c(plot_Ds, rev(plot_Ds)),
c(pad_inner_means + 2 * pad_inner_ses, rev(pad_inner_means - 2 * pad_inner_ses)),
col = c_light, border = NA)
lines(plot_Ds, pad_inner_means, col=c_light_highlight, lwd=2)
polygon(c(plot_Ds, rev(plot_Ds)),
c(pad_outer_means + 2 * pad_outer_ses, rev(pad_outer_means - 2 * pad_outer_ses)),
col = c_dark, border = NA)
lines(plot_Ds, pad_outer_means, col=c_dark_highlight, lwd=2)

# Plot ratio of probabilities verses dimension
plot(1, type="n", main="",
xlim=c(head(Ds, 1) - 0.5, tail(Ds, 1) + 0.5), xlab="Dimension",
ylim=c(0, 10), ylab="Ratio of Outer verses Inner Probability")

lines(plot_Ds, pad_outer_means / pad_inner_means, col=c_dark_highlight, lwd=2)

############################################################
#
# How does the distance between two sampled points behave
# as the dimensionality of the box increases?
#
############################################################

delta_means <- 0 * Ds
delta_ses <- 0 * Ds
delta_quantiles <- matrix(data = 0, nrow = tail(Ds, 1), ncol = 9)

for (D in Ds) {
  # Distances between two sampled points
  delta_samples <- rep(0, N)
  
  for (n in 1:N) {
    # Sample two points
    x1 <- runif(D, -3, 3)
    x2 <- runif(D, -3, 3)
    
    # Compute distance between them
    delta_samples[n] <- sqrt(sum( (x1 - x2)**2 ))
  }
  
  # Estimate average distance
  s <- compute_mc_stats(delta_samples)
  delta_means[D] <- s[1]
  delta_ses[D] <- s[2]
  
  # Estimate distance quantiles
  delta_quantiles[D,] <- quantile(delta_samples, probs=quant_probs)
}

# Plot average distance between points verses dimension
pad_delta_means <- do.call(cbind, lapply(idx, function(n) delta_means[n]))
pad_delta_ses <- do.call(cbind, lapply(idx, function(n) delta_ses[n]))

plot(1, type="n", main="",
xlim=c(head(Ds, 1) - 0.5, tail(Ds, 1) + 0.5), xlab="Dimension",
ylim=c(0, 10), ylab="Average Distance Between Points")

polygon(c(plot_Ds, rev(plot_Ds)),
c(pad_delta_means + 2 * pad_delta_ses, rev(pad_delta_means - 2 * pad_delta_ses)),
col = c_light, border = NA)
lines(plot_Ds, pad_delta_means, col=c_dark, lwd=2)

# Plot distance quantiles verses dimension
pad_delta_quantiles <- do.call(cbind, lapply(idx, function(n) delta_quantiles[n, 1:9]))

plot(1, type="n", main="",
xlim=c(head(Ds, 1) - 0.5, tail(Ds, 1) + 0.5), xlab="Dimension",
ylim=c(0, 10), ylab="Distance Between Points")

polygon(c(plot_Ds, rev(plot_Ds)), c(pad_delta_quantiles[1,], rev(pad_delta_quantiles[9,])),
col = c_light, border = NA)
polygon(c(plot_Ds, rev(plot_Ds)), c(pad_delta_quantiles[2,], rev(pad_delta_quantiles[8,])),
col = c_light_highlight, border = NA)
polygon(c(plot_Ds, rev(plot_Ds)), c(pad_delta_quantiles[3,], rev(pad_delta_quantiles[7,])),
col = c_mid, border = NA)
polygon(c(plot_Ds, rev(plot_Ds)), c(pad_delta_quantiles[4,], rev(pad_delta_quantiles[6,])),
col = c_mid_highlight, border = NA)
lines(plot_Ds, pad_delta_quantiles[5,], col=c_dark, lwd=2)

############################################################
#
# How does the distance from a Gaussian sample and the
# Gaussian mode behave as the dimensionality increases?
#
############################################################

r_means <- 0 * Ds
r_ses <- 0* Ds
r_quantiles <- matrix(data = 0, nrow = tail(Ds, 1), ncol = 9)

for (D in Ds) {
  # Distance from Gaussian samples to mode at zero
  r_samples <- rep(0, N)
  
  for (n in 1:N) {
    # Sample point
    x <- rnorm(D, 0, 1)
    
    # Compute distance from point to mode at zero
    r_samples[n] <- sqrt(sum(x**2))
  }
  
  # Estimate average distance
  s <- compute_mc_stats(r_samples)
  r_means[D] <- s[1]
  r_ses[D] <- s[2]
  
  # Estimate distance quantiles
  r_quantiles[D,] <- quantile(r_samples, probs=quant_probs)
}

# Plot average distance from mode verses dimension
pad_r_means <- do.call(cbind, lapply(idx, function(n) r_means[n]))
pad_r_ses <- do.call(cbind, lapply(idx, function(n) r_ses[n]))

plot(1, type="n", main="",
xlim=c(head(Ds, 1) - 0.5, tail(Ds, 1) + 0.5), xlab="Dimension",
ylim=c(0, 4), ylab="Average Distance From Mode")

polygon(c(plot_Ds, rev(plot_Ds)),
c(pad_r_means + 2 * pad_r_ses, rev(pad_r_means - 2 * pad_r_ses)),
col = c_light, border = NA)
lines(plot_Ds, pad_r_means, col=c_dark, lwd=2)

# Plot distance quantiles verses dimension
pad_r_quantiles <- do.call(cbind, lapply(idx, function(n) r_quantiles[n, 1:9]))

plot(1, type="n", main="",
xlim=c(head(Ds, 1) - 0.5, tail(Ds, 1) + 0.5), xlab="Dimension",
ylim=c(0, 4), ylab="Distance From Mode")

polygon(c(plot_Ds, rev(plot_Ds)), c(pad_r_quantiles[1,], rev(pad_r_quantiles[9,])),
col = c_light, border = NA)
polygon(c(plot_Ds, rev(plot_Ds)), c(pad_r_quantiles[2,], rev(pad_r_quantiles[8,])),
col = c_light_highlight, border = NA)
polygon(c(plot_Ds, rev(plot_Ds)), c(pad_r_quantiles[3,], rev(pad_r_quantiles[7,])),
col = c_mid, border = NA)
polygon(c(plot_Ds, rev(plot_Ds)), c(pad_r_quantiles[4,], rev(pad_r_quantiles[6,])),
col = c_mid_highlight, border = NA)
lines(plot_Ds, pad_r_quantiles[5,], col=c_dark, lwd=2)

# Plot residual quantiles verses dimension
plot(1, type="n", main="",
xlim=c(head(Ds, 1) - 0.5, tail(Ds, 1) + 0.5), xlab="Dimension",
ylim=c(-1, 1), ylab="Residual from Median Distance")

polygon(c(plot_Ds, rev(plot_Ds)),
c(pad_r_quantiles[1,] - pad_r_quantiles[5,], rev(pad_r_quantiles[9,] - pad_r_quantiles[5,])),
col = c_light, border = NA)
polygon(c(plot_Ds, rev(plot_Ds)),
c(pad_r_quantiles[2,] - pad_r_quantiles[5,], rev(pad_r_quantiles[8,] - pad_r_quantiles[5,])),
col = c_light_highlight, border = NA)
polygon(c(plot_Ds, rev(plot_Ds)),
c(pad_r_quantiles[3,] - pad_r_quantiles[5,], rev(pad_r_quantiles[7,] - pad_r_quantiles[5,])),
col = c_mid, border = NA)
polygon(c(plot_Ds, rev(plot_Ds)),
c(pad_r_quantiles[4,] - pad_r_quantiles[5,], rev(pad_r_quantiles[6,] - pad_r_quantiles[5,])),
col = c_mid_highlight, border = NA)
lines(plot_Ds, pad_r_quantiles[5,] - pad_r_quantiles[5,], col=c_dark, lwd=2)

############################################################
#
# What is the probability of a Gaussian sample falling
# into a spherical neighborhood around the mode at zero?
#
############################################################

prob_means <- 0 * Ds
prob_ses <- 0 * Ds

R <- 1

for (D in Ds) {
  # Does the sample fall into the spherical neighborhood?
  is_central_samples <- rep(0, N)
  
  for (n in 1:N) {
    # Sample a new point
    x <- rnorm(D, 0, 1)
    
    # Compute radial distance from mode
    r <- sqrt(sum(x**2))
    
    # Check if sample is contained within spherical neighborhood
    if (r < R)
    is_central_samples[n] <- 1
  }
  
  # Estimate probability of falling into spherical neighborhood
  s <- compute_mc_stats(is_central_samples)
  prob_means[D] <- s[1]
  prob_ses[D] <- s[2]
}

# Plot inclusion probability verses dimension
pad_means <- do.call(cbind, lapply(idx, function(n) prob_means[n]))
pad_ses <- do.call(cbind, lapply(idx, function(n) prob_ses[n]))

plot(1, type="n", main="",
xlim=c(head(Ds, 1) - 0.5, tail(Ds, 1) + 0.5), xlab="Dimension",
ylim=c(0, 0.7), ylab="Inclusion Probability")

polygon(c(plot_Ds, rev(plot_Ds)), c(pad_means + 2 * pad_ses, rev(pad_means - 2 * pad_ses)),
col = c_light, border = NA)
lines(plot_Ds, pad_means, col=c_dark, lwd=2)
