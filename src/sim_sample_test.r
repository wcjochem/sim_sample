
# 
# Simulate spatially-correlated density maps
# sim_sample_test.r
#
# Demonstration program for simulating spatially-correlated data
# and testing a range of sampling techniques to estimate population.
#
# Author: Chris Jochem (W.C.Jochem@soton.ac.uk)
# March 2018; update: 04 May 2018
#
# Random field simulation techniques based on work by:
# Jérôme Guélat (https://rpubs.com/jguelat/autocorr2)
# Dorazio and Karanth (2017). "A hierarchical model for estimating 
#    the spatial distribution and abundance of animals detected 
#    by continuous-time recorders." PLoS ONE, 12(5): e0176966.
#

library(raster)
library(gstat)

# simulation parameters
nsims <- 10 # number of spatial fields to repeat
regionsize <- 50 # dim of square region
phi <- 0.05 # smoothness parameter

# abundance mdoel parameters
beta0 <- .4 # intercept
beta1 <- 3 # elevation
beta2 <- -1.5 # quadratic
beta3 <- -1 # trend


# Function to draw samples from a multivariate normal distribution
rmvn <- function(n, mu = 0, V = matrix(1)) {
  p <- length(mu)
  if (any(is.na(match(dim(V), p)))) 
    stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(n * p), ncol = p) %*% D + rep(mu, rep(n, p)))
}

## create covariate layers
# elevation
data("volcano") # R dataset for example
# convert matrix to raster
# clipping original down to size
elev <- raster(volcano[10:59,5:54],
               xmn=0, xmx=regionsize, ymn=0, ymx=regionsize)
elev <- scale(elev) # z-scores
  plot(elev)
  
# trend surface - smoothly varying over space
trend <- raster(elev) # copy elevation raster extent
# get cells to calculate surface values
locns <- xyFromCell(trend, 1:ncell(trend))
# calculate smoothly varying value by x,y coords
vals <- -0.01*locns[,'x'] - 0.01*locns[,'y'] + 0.01*locns[,'x']*locns[,'y']
# update raster values
values(trend) <- vals
trend <- scale(trend) # z-scores
  plot(trend)

## create spatial fields
# pre-allocate storage
countfields <- vector("list", length=nsims)
# loop to create abundance datasets
# NOTE this could be slow for large nsims and/or regionsize
for(i in 1:nsims){
  print(i)
  set.seed(i)
  # Set up a square study region
  simgrid <- expand.grid(1:regionsize, 1:regionsize)
  n <- nrow(simgrid)
  # Set up distance matrix
  distance <- as.matrix(dist(simgrid))
  # Generate random variable
  X <- rmvn(1, rep(0, n), exp(-phi * distance))
  
  spfield <- rasterFromXYZ(cbind(simgrid[, 1:2] - 0.5, X))
    # plot(spfield)
  # abundance model
  # lambda <- exp(beta0 + beta1*values(elev) + beta2*values(trend) + beta3*(values(elev) * values(trend)) + values(spfield))
  # lambda <- exp(beta0 + beta1*values(elev) + beta2*values(trend) + values(spfield))
  lambda <- exp(beta0 + beta1*values(elev) + beta2*values(elev)^2 + beta3*values(trend) + values(spfield))
    # plot(values(elev), lambda, cex=0.5)
  # generate counts from the simulated intensity
  counts <- rpois(n, lambda)
    sum(counts)
  # gridded abundance counts
  counts <- rasterFromXYZ(cbind(coordinates(elev), counts))
    # plot(counts)
  # store the simulations
  countfields[[i]] <- counts
}
# can save/output the list of simulated fields for replication study

# visualise the results
plot(counts[[1]]) # first example

