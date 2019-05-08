# 
# Analysis code for simulated population surfaces
# sim_analysis.r
#
#
# Author: Chris Jochem (W.C.Jochem@soton.ac.uk)
#

library(INLA)
library(Metrics)


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

# create spatial mesh for geostatistical models
fake.locations <- matrix(c(0,0,r_width,r_height,0,r_height,r_width,0), nrow=4, byrow=T)
sim_mesh <- inla.mesh.2d(loc=fake.locations, max.edge=c(0.4, 2)) # matches mesh used to generate data
  
# sampling  
nsamples <- 100 # number of draws of each size and type
sampleszs <- c(50, 100, 150, 200) # range of sample sizes

# get list of simulated surfaces
sim_files <- list.files("./out_sim", pattern=glob2rx("sim_surface_*.rds"), recursive=FALSE, full.names=TRUE)
# set up storage for results
out_length <- length(sim_files) * 50 * nsamples * 5 * length(sampleszs) * 4 # files x fields x draws x sample tech x sizes x models

all_outputs <- data.frame(field=vector("numeric", length=out_length),
                          surface=vector("numeric", length=out_length),
                          sample=vector("numeric", length=out_length),
                          samp_type=vector("character", length=out_length),
                          size=vector("numeric", length=out_length),
                          model=vector("character", length=out_length),
                          rmse_a=vector("numeric", length=out_length), # _a are within sample area predictions
                          mape_a=vector("numeric", length=out_length),
                          pctcov_a=vector("numeric", length=out_length),
                          pr2_a=vector("numeric", length=out_length),
                          rmse_b=vector("numeric", length=out_length), # _b are out-of-sample predictions
                          mape_b=vector("numeric", length=out_length),
                          pctcov_b=vector("numeric", length=out_length),
                          pr2_b=vector("numeric", length=out_length),
                          stringsAsFactors=F
                         )


## Main processing loop ##
r <- 1 # index to store the results
# loop and load all files
for(f in sim_files){
  print(f)
  sim_surfaces <- readRDS(f) # load the simulated population and covariates
  nsurfaces <- length(sim_surfaces) # number of simulated populations
  truerange <- as.numeric(gsub('^.*e_\\s*|\\s*.rds.*$', '', f)) # get range of spatial field
  
  # loop each simulated surface
  for(s in 1:nsurfaces){
    # get simulated layers
    truepop <- sim_surfaces[[s]]$pop
    sett <- sim_surfaces[[s]]$sett
    cov <- sim_surfaces[[s]]$cov
    
    # create the domain dataset
    domain <- as.data.frame(truepop, xy=TRUE)
    names(domain)[3] <- "pop"
    # make cell numbers for linking to other grids
    domain$cnum <- cellFromXY(truepop, domain[,c("x","y")])
    # extract other values
    domain$sett <- extract(sett, domain[,c("x","y")])
    domain$cov <- extract(cov, domain[,c("x","y")])
    # get region information
    domain$adm1 <- extract(strat1, domain[,c("x","y")])
    domain$adm2 <- extract(strat2, domain[,c("x","y")])
    # exclude areas with no population
    domain <- domain[domain$pop > 0,]
    
    # split domain into obs/pred
    obs <- domain[domain$x<10,]
    pred <- domain[domain$x>=10,]
    # get aggregate population totals
    totobspop <- sum(obs$pop)
    totpredpop <- sum(pred$pop)
    
    # loop sample size
    for(sz in sampleszs){
      # repeated samples of different types
      for(sp in 1:nsamples){
        # simple random sample
        srs <- sample(1:nrow(obs), size=sz, replace=F)
        # extract values
        srs <- obs[srs,]
        
        # simple mean population per pixel
        obs$mean <- mean(srs$pop)
        pred$mean <- mean(srs$pop)
        # store results
        all_outputs[r, c("field","surface","sample","size")] <- c(truerange, s, sp, sz)
        all_outputs[r, c("samp_type", "model")] <- c("srs","mean")
        all_outputs[r, c("rmse_a","mape_a","pctcov_a","pr2_a")] <- c(rmse(obs$pop, obs$mean), mape(obs$pop, obs$mean), NA, cor(obs$pop, obs$mean)^2)
        all_outputs[r, c("rmse_b","mape_b","pctcov_b","pr2_b")] <- c(rmse(pred$pop, pred$mean), mape(pred$pop, pred$mean), NA, cor(pred$pop, pred$mean)^2)
        r <- r + 1 # increment counter
        
        # geostatistical model
        
        
      } # end multiple sample draws loop
      
    } # end sample sizes loop
    
  } # end loop of simulated surfaces
  
} # end loop of files
