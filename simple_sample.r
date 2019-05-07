# 
# Simulate spatially-correlated density maps
# sim_sample.r
#
# Demonstration program for simulating spatially-correlated data
# and testing a range of sampling techniques to estimate population.
#
# Author: Chris Jochem (W.C.Jochem@soton.ac.uk)
# Simulation set-up based on: https://haakonbakka.bitbucket.io/btopic122.html
#

library(INLA)
library(fields)
library(ggplot2)
library(viridisLite)
library(raster)
library(randomcoloR)


rm(list=ls())
set.seed(1126)
# inla.seed <- 1126

####
# helper function to plot spatial field
# source: https://haakonbakka.bitbucket.io/btopic108.html
local.plot.field <- function(field, mesh, xlim=c(0,10), ylim=c(0,10), ...){
  stopifnot(length(field) == mesh$n)
  
  proj <- inla.mesh.projector(mesh, xlim=xlim, ylim=ylim, dims=c(500, 500))
  field.proj <- inla.mesh.project(proj, field)
  
  n.col <- 20
  image.plot(list(x=proj$x, y=proj$y, z=field.proj), xlim=xlim, ylim=ylim, col=plasma(n.col), nlevel=n.col+1, ...)
}

# load code to create settlement type grids
source("./make_settlement.r")


## simulations set-up ##
# dimensions
r_width <- 20
r_height <- 20
# create spatial units 
# make 5 horizontal strata (naive stratification for comparison)
# see: https://gis.stackexchange.com/questions/34895/create-zonal-grid-in-r
rows <- 50; cols <- 200; n <- 200
strat1 <- outer(1:n, 1:n, 
               function(i,j) (i-1) %/% rows * ((n+1) %/% cols) + (j-1) %/% cols + 1)
strat1 <- raster(strat1, xmn=0, xmx=20, ymn=0, ymx=20)
# nstrat <- length(unique(values(strat)))
  plot(strat1, col=randomColor(4), legend=F)
  abline(v=10)
# level2 local areas
rows <- 25; cols <- 25; n <- 200
strat2 <- outer(1:n, 1:n, 
               function(i,j) (i-1) %/% rows * ((n+1) %/% cols) + (j-1) %/% cols + 1)
strat2 <- raster(strat2, xmn=0, xmx=20, ymn=0, ymx=20)
  plot(strat2, col=randomColor(64), legend=F)
  abline(v=10)
  
## population model set-up ##
# true values
# sigma.u <- 1.5
# range <- 3 # 1, 3, 5
# 
# # make mesh
# fake.locations <- matrix(c(0,0,r_width,r_height,0,r_height,r_width,0), nrow=4, byrow=T)
# # mesh.sim <- inla.mesh.2d(loc=fake.locations, max.edge=c(0.4, 2), offset=c(1,6))
# mesh.sim <- inla.mesh.2d(loc=fake.locations, max.edge=c(0.4, 2))
# 
#   plot(mesh.sim)
#   points(fake.locations, pch=16)
#   axis(1); axis(2)
#   
# # simulate field  
# spde <- inla.spde2.pcmatern(mesh.sim, prior.range=c(.5, .5), prior.sigma=c(.5, .5))
# 
# Qu <- inla.spde.precision(spde, theta=c(log(range), log(sigma.u)))
# u <- inla.qsample(n=1, Q=Qu, seed=inla.seed)
# u <- u[ ,1]
# 
#   local.plot.field(u, mesh.sim, xlim=c(0,r_width), ylim=c(0,r_height))
#   len <- range
#   arrows(5-0.5*len, .5, 5+0.5*len, .5, length=0.05, angle=90, code=3, lwd=3)
#   
# # simulate data at locations
# loc.data <- as.matrix(expand.grid(seq(0, r_width, .1), seq(0, r_height, .1))) + .1
# n <- nrow(loc.data)
# # create a raster covering the study area
# r <- raster(nrows=r_height, ncols=r_width, xmn=0, xmx=r_width, ymn=0, ymx=r_height, resolution=.1, crs=NULL)
# loc.data <- coordinates(r)
# n <- nrow(loc.data)
# 
# # project spatial field to locations
# A <- inla.spde.make.A(mesh=mesh.sim, loc=loc.data)
# u <- drop(A %*% u)
# 
# # create a covariate
# x <- runif(n)-0.5 # centred at zero
# cov_r <- r # create a blank raster
# values(cov_r) <- x # write in values
#   plot(cov_r)
# 
# # construct linear predictor
# beta <- c(1, 1.5) # true coefficients 
# 
# lin.pred <- beta[1] + beta[2]*x + u
# # observation
# y <- rpois(n, exp(lin.pred))
# 
# y_rast <- r # store pop as raster
# values(y_rast) <- y
#   plot(y_rast) # example gridded population
# # calculate totals per area (test)
# zonal(y_rast, strat2)
# zonal(y_rast, strat1)
# sum(y)
# 
# # derive a settlement "covariate" based on population centers/density
# sett <- make_settlement(y_rast)
#   plot(sett, col=randomColor(6)) # raster plotting
#   abline(v=10)
# 
# zonal(y_rast, sett)  


## loop to build synthetic populations ##
# set up values
ranges <-  c(1.5,3,5) # varying spatial field
#range <- 3
sigma.u <- 1
nsims <- 5 # number of population per spatial field

fake.locations <- matrix(c(0,0,r_width,r_height,0,r_height,r_width,0), nrow=4, byrow=T)
mesh.sim <- inla.mesh.2d(loc=fake.locations, max.edge=c(0.4, 2))
# simulate field  
spde <- inla.spde2.pcmatern(mesh.sim, prior.range=c(.5, .5), prior.sigma=c(.5, .5))
# simulate data at locations
loc.data <- as.matrix(expand.grid(seq(0, r_width, .1), seq(0, r_height, .1))) + .1
n <- nrow(loc.data)
# create a raster covering the study area
r <- raster(nrows=r_height, ncols=r_width, xmn=0, xmx=r_width, ymn=0, ymx=r_height, resolution=.1, crs=NULL)
loc.data <- coordinates(r)
n <- nrow(loc.data)


# storage
simdata <- vector(mode="list", length=nsims)

# looping
for (range in ranges){ # loop and vary spatial range
  print(paste0("Range: ", range))
  print(" Simulations")

  for (s in 1:nsims){ # loop number of simulated pops
    print(paste0(" ", s))
    inla.seed <- s # fix seed per simulation
    # simulate spatial field
    Qu <- inla.spde.precision(spde, theta=c(log(range), log(sigma.u)))
    u <- inla.qsample(n=1, Q=Qu, seed=inla.seed)
    u <- u[ ,1]
      
    # simulate data at locationss
    loc.data <- as.matrix(expand.grid(seq(0, r_width, .1), seq(0, r_height, .1))) + .1
    n <- nrow(loc.data)
    # create a raster covering the study area
    r <- raster(nrows=r_height, ncols=r_width, xmn=0, xmx=r_width, ymn=0, ymx=r_height, resolution=.1, crs=NULL)
    loc.data <- coordinates(r)
    n <- nrow(loc.data)
    
    # project spatial field to locations
    A <- inla.spde.make.A(mesh=mesh.sim, loc=loc.data)
    u <- drop(A %*% u)
    
    # create a covariate
    x <- runif(n)-0.5 # centred at zero
    cov_r <- r # create a blank raster
    values(cov_r) <- x # write in values
      # plot(cov_r)
    
    # construct linear predictor
    beta <- c(1, 1.5) # true coefficients 
    
    lin.pred <- beta[1] + beta[2]*x + u
    # observation
    y <- rpois(n, exp(lin.pred))
    
    y_rast <- r # store pop as raster
    values(y_rast) <- y
    
    # derive settlement "types" from population
    sett <- make_settlement(y_rast)
  
    simdata[[s]] <- list("pop"=y_rast, "cov"=cov_r, "sett"=sett)
  } # end pop simulation loop
  # save simulation outputs
  saveRDS(simdata, file=paste0("./out_sim/sim_surface_", range,".rds"))
} # end loop over range values








