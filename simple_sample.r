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
inla.seed <- 1126

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

# true values
sigma.u <- 1.5
range <- 3 # 1, 3, 5

# make mesh
fake.locations <- matrix(c(0,0,r_width,r_height,0,r_height,r_width,0), nrow=4, byrow=T)
# mesh.sim <- inla.mesh.2d(loc=fake.locations, max.edge=c(0.4, 2), offset=c(1,6))
mesh.sim <- inla.mesh.2d(loc=fake.locations, max.edge=c(0.4, 2))

  plot(mesh.sim)
  points(fake.locations, pch=16)
  axis(1); axis(2)
  
# simulate field  
spde <- inla.spde2.pcmatern(mesh.sim, prior.range=c(.5, .5), prior.sigma=c(.5, .5))

Qu <- inla.spde.precision(spde, theta=c(log(range), log(sigma.u)))
u <- inla.qsample(n=1, Q=Qu, seed=inla.seed)
u <- u[ ,1]

  local.plot.field(u, mesh.sim, xlim=c(0,r_width), ylim=c(0,r_height))
  len <- range
  arrows(5-0.5*len, .5, 5+0.5*len, .5, length=0.05, angle=90, code=3, lwd=3)
  
# simulate data at locations
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
  plot(cov_r)

# construct linear predictor
beta = c(1, 1.5) # true coefficients 

lin.pred <- beta[1] + beta[2]*x + u
# observation
y <- rpois(n, exp(lin.pred))

y_rast <- r # store pop as raster
values(y_rast) <- y
  plot(y_rast)
  
