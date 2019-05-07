# 
# Analysis code for simulated population surfaces
# sim_analysis.r
#
#
# Author: Chris Jochem (W.C.Jochem@soton.ac.uk)
#

library(INLA)


rm(list=ls())
set.seed(1126)

# set-up
# construct study area and strata
r_width <- 20
r_height <- 20
# create spatial units 
# make 5 horizontal strata (naive stratification for comparison)
# see: https://gis.stackexchange.com/questions/34895/create-zonal-grid-in-r
rows <- 50; cols <- 200; n <- 200
strat1 <- outer(1:n, 1:n, 
               function(i,j) (i-1) %/% rows * ((n+1) %/% cols) + (j-1) %/% cols + 1)
strat1 <- raster(strat1, xmn=0, xmx=20, ymn=0, ymx=20)
# level2 local areas
rows <- 25; cols <- 25; n <- 200
strat2 <- outer(1:n, 1:n, 
               function(i,j) (i-1) %/% rows * ((n+1) %/% cols) + (j-1) %/% cols + 1)
strat2 <- raster(strat2, xmn=0, xmx=20, ymn=0, ymx=20)
  
# sampling  
nsamples <- 100 # number of draws of each size and type
sampleszs <- c(50, 100, 150, 200) # range of sample sizes

# get list of simulated surfaces
sim_files <- list.files("./out_sim", pattern=glob2rx("sim_surface_*.rds"), recursive=FALSE, full.names=TRUE)
# loop and load all files
for(f in sim_files){
  print(f)
  sim_surfaces <- readRDS(f) # load the simulated population and covariates
  nsurfaces <- length(sim_surfaces) # number of simulated populations
  
  # loop each simulated surface
  for(s in 1:nsurfaces){
    # get simulated layers
    truepop <- sim_surfaces[[s]]$pop
    sett <- sim_surfaces[[s]]$sett
    cov <- sim_surfaces[[s]]$cov
    
  } # end loop of simulated surfaces
  
} # end loop of files
