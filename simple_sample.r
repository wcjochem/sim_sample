# 
# Simulate spatially-correlated density maps
# sim_sample.r
#
# Demonstration program for simulating spatially-correlated data
# and testing a range of sampling techniques to estimate population.
#
# Author: Chris Jochem (W.C.Jochem@soton.ac.uk)
#
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

