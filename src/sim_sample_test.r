
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
library(ggplot2)
library(rasterVis)
library(gridExtra)
library(splitstackshape)

# simulation parameters
nsims <- 10 # number of spatial fields to repeat
regionsize <- 50 # dim of square region
phi <- 0.05 # smoothness parameter
ndraws <- 100 # number of repeated samples

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
plot(countfields[[1]]) # first example
# create list of all counts
count_plots <- lapply(1:nsims, function(x){
  gplot(countfields[[x]]) +
    geom_tile(aes(fill= value)) + 
    scale_fill_gradient(low="white", high="red") + 
    theme_bw() +
    coord_equal()
})
  # count_plots
# arrange all plots
marrangeGrob(count_plots, nrow=5, ncol=2)


## implement sampling strategies
# different sample sizes
sampsz <- c(25,50,75,100,125,150,200,250,300,400,500) # CHANGE HERE

# preallocate storage for the output
pred_errs <- matrix(data=cbind(rep(1:nsims, each=length(sampsz)*ndraws),
                               rep(sampsz, each=ndraws),
                               NA, NA), # create multiple cols for error metrics
                    nrow=nsims*length(sampsz)*ndraws, ncol=4)
# Loop and evaluate all
r <- 1 # counter to fill in the output
for(i in 1:nsims){ # loop over different simulated populations
  print(i)
  
  count <- countfields[[i]] # simulated pop
  # total "true" population
  totpop <- sum(values(count))
  
  # create binary 'settled' area
  settle <- count # copy raster
  settle[settle>0] <- 1
  # total number of settled cells (need later for proportional allocation)
  totsettle <- sum(values(settle))
  # create a domain to select sampling locations
  # sampling from within settlement (assumes perfect residential classification)
  # extract records to a data.frame (ok for small simulations)
  domain <- as.data.frame(settle, xy=TRUE)
  domain <- domain[domain$counts>0,]
  # extract observed (simulated) population values for all locations
  domain$counts <- extract(count, domain[,c("x","y")])
  # "true" average population per pixel
    mean(domain$counts)
    
  # create a high pop sub-domain area to estimate totals
  # analogous to having a capital city within a regional sampling domain
  hicell <- cellFromXY(count, domain[which.max(domain$counts),c("x","y")])
  row <- rowFromCell(count, hicell) # get row/column of the high population area
  col <- colFromCell(count, hicell)

  csz <- 4 # set sub-area size (square: csz+1 x csz+1 cells)
  stopifnot(regionsize >= 2*csz + 1) # double-check no small regions
  rcols <- lcols <- trows <- brows <- csz # default square sub-region
  # logic checks to prevent subdomain extending beyond border
  if(col+csz > regionsize){
    rcols <- regionsize - col
    lcols <- lcols + (csz - rcols)
  }
  if(col-csz <= 0){
    lcols <- col - 1
    rcols <- csz + (csz-lcols)
  }
  # repeat checks for rows
  if(row+csz > regionsize){
    brows <- regionsize - row
    trows <- csz + (csz - brows)
  }
  if(row-csz <= 0){
    trows <- row - 1
    brows <- csz + (csz - trows)
  }
  # create sub-region extent
  subreg <- extent(count, row-trows, row+brows, col-lcols, col+rcols)
    plot(count); points(domain[which.max(domain$counts),c("x","y")]); plot(subreg, add=T) # example plot
  # extract subregion to new population grid
  subpop <- crop(count, subreg)
    # sum(values(subpop)) # quick count of total pop
  
  ## stratify the settled area - multiple methods
  # make 5 horizontal strata (naive stratification for comparison)
  # see: https://gis.stackexchange.com/questions/34895/create-zonal-grid-in-r
  rows <- 10; cols <- 50; n <- 50
  strat <- outer(1:n, 1:n, 
                 function(i,j) (i-1) %/% rows * ((n+1) %/% cols) + (j-1) %/% cols + 1)
  strat <- raster(strat, xmn=0, xmx=50, ymn=0, ymx=50)
  nstrat <- length(unique(values(strat)))
    # plot(count); plot(strat, add=T, alpha=.5) # example
  # create strata ID for sampling domain
  domain$simp_strat <- as.character(extract(strat, domain[,c("x","y")]))
  
  # population-based strata
  # create weights from the aggregated, average pop counts
  count_agg5 <- aggregate(count, fact=5, fun=mean)
    # plot(count_agg5)
  # extract values to the sample domain
  domain$pop_wgt <- extract(count_agg5, domain[,c("x","y")])
  
  for(sz in sampsz){ # vary total sample size
    print(sz)
  # repeated sample realisations
    for(d in 1:ndraws){ 
      # set.seed(d) # ?
    # sampling methods:
    ## simple random sample ##
      srs <- sample(1:nrow(domain), size=sz, replace=F)
      # extract values of sampled points
      srs <- domain[srs,]
      # sample mean pop per pixel
      samp_mean <- mean(srs$counts)
      
    ## stratified random sample - equal weight ##
      # equal sample size per stratum
      szs <- round(sz / nstrat)
      # draw sample
      strs <- data.frame(stratified(domain, "simp_strat", szs))
      # TO-DO: calculate stratified mean values of [counts]
      
   ## stratified random sample - prop to settled area ##
      wgt <- zonal(settle, strat, sum)
      wgt <- wgt[,'value'] / totsettle
      # sample sizes - note rounding may causes some different totals
      szs <- setNames(round(wgt*sz), unique(values(strat)))
      # draw sample
      strs <- data.frame(stratified(domain, "simp_strat", szs))
      # TO-DO: calculate stratified mean values
      
   ## sample weighted by approximate population density
      srs_pwgt <- sample(1:nrow(domain), size=sz, replace=F, prob=domain$pop_wgt)
      # extract values from sampled points
      srs_pwgt <- domain[srs_pwgt,]
      # sample mean pop per pixel
      
      
      # TO-DO: calculate and store the comparison of total population
      # TO-DO: compare prediction for sub-region population
      # pred_errs[r,3] <- ERROR calculation storage
      # pred_errs[r, 4] <- ERROR calculation storage
      r <- r+1
    }
  }
}

# process results
