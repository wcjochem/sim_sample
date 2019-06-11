# 
# Analysis code for 1 simulated population
# sim_analysis_1.r
#
#
# Author: Chris Jochem (W.C.Jochem@soton.ac.uk)
#

library(raster)
library(INLA)
library(Metrics)
library(BalancedSampling)
library(splitstackshape)


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
sim_mesh <- inla.mesh.2d(loc=fake.locations, max.edge=c(1, 4)) 
  plot(sim_mesh, main="")

# load simulations
# get list of simulated surfaces
sim_files <- list.files("./out_sim", pattern=glob2rx("sim_surface_*.rds"), recursive=FALSE, full.names=TRUE)
f <- sim_files[[2]]
 print(f)
sim_surfaces <- readRDS(f) # load the simulated population and covariates
truerange <- as.numeric(gsub('^.*e_\\s*|\\s*.rds.*$', '', f)) # get range of spatial field
# get simulated layers
truepop <- sim_surfaces[[3]]$pop
sett <- sim_surfaces[[3]]$sett
cov <- sim_surfaces[[3]]$cov

  par(mfrow=c(1,3))
  plot(truepop)
  plot(sett)
  plot(cov)
  par(mfrow=c(1,1))
  
# find study sites ("towns") and create extents
  plot(SpatialPoints(matrix(c(16.5, 11), ncol=2, byrow=T)), add=T)
csz <- 9 # area size
pt1 <- cellFromXY(truepop, c(7.5, 7.5))
loc1r <- rowFromCell(truepop, pt1)
loc1c <- colFromCell(truepop, pt1)
loc1x <- extent(truepop, loc1r-csz, loc1r+csz, loc1c-csz, loc1c+csz)

pt2 <- cellFromXY(truepop, c(16, 15.5))
loc2r <- rowFromCell(truepop, pt2)
loc2c <- colFromCell(truepop, pt2)
loc2x <- extent(truepop, loc2r-csz, loc2r+csz, loc2c-csz, loc2c+csz)

pt3 <- cellFromXY(truepop, c(16.5, 11))
loc3r <- rowFromCell(truepop, pt3)
loc3c <- colFromCell(truepop, pt3)
loc3x <- extent(truepop, loc3r-csz, loc3r+csz, loc3c-csz, loc3c+csz)

pt4 <- cellFromXY(truepop, c(4, 5))
loc4r <- rowFromCell(truepop, pt4)
loc4c <- colFromCell(truepop, pt4)
loc4x <- extent(truepop, loc4r-csz, loc4r+csz, loc4c-csz, loc4c+csz)


  plot(truepop, asp=1)
  plot(sett, asp=1)
  abline(v=10, lty=2)
  plot(loc1x, add=T, lwd=2, col="blue")
  plot(loc2x, add=T, lwd=2, col="red")
  plot(loc3x, add=T, lwd=2, col="orange")
  plot(loc4x, add=T, lwd=2, col="purple")

# create datasets for modelling
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
# get sub regions
domain[domain$cnum %in% cellsFromExtent(truepop, loc1x),"town"] <- 1
domain[domain$cnum %in% cellsFromExtent(truepop, loc2x),"town"] <- 2
domain[domain$cnum %in% cellsFromExtent(truepop, loc3x),"town"] <- 3
domain[domain$cnum %in% cellsFromExtent(truepop, loc4x),"town"] <- 4
# exclude areas with no population
domain <- domain[domain$pop > 0,]
  by(domain$pop, domain$town, sum)

# split domain into obs/pred
obs <- domain[domain$x<10,]
# create unique row ID for sampling
obs$rid <- 1:nrow(obs)
pred <- domain[domain$x>=10,]


## samples ##
sz <- 50 # fixed sample size
# simple random sample
srs <- sample(1:nrow(obs), size=sz, replace=F)
# population weighted
wts <- (obs$pop/sum(obs$pop))*sz
pwgt <- sample(1:nrow(obs), size=sz, prob=wts, replace=F)
# stratified random by type - proportional to area

# systematic random by admin2
szs <- trunc(sz / length(unique(obs$adm2)))
sys <- data.frame(stratified(obs, "adm2", szs))$rid
# redistribute other samples randomly
rem <- sz %% length(unique(obs$adm2))
# get blocks
extra_adm <- sample(unique(obs$adm2), size=rem, replace=F)
# sample a pixel within admin unit
extra <- data.frame(stratified(obs[-sys,], "adm2", 1, select=list(adm2=extra_adm)))$rid
sys <- c(sys,extra) # combine

# spatially-balanced
# population-weighted
bs <- lcube(prob=wts, 
            Xspread=as.matrix(obs[,c("x","y")]),
            Xbal=as.matrix(cbind(wts,obs[,c("cov","sett","x","y")]))
           )

# set samples
obs[srs,"srs"] <- TRUE
obs[pwgt, "pwt"] <- TRUE
obs[bs,"bal"] <- TRUE

sampdf <- obs[samps,]
pred_a <- obs[-samps,] # in-domain predictions











# balanced spatial sample
# bs <- lcubestratified(rep(50/nrow(obs),nrow(obs)), 
# bs <- lcubestratified(rep(50/nrow(obs),nrow(obs)), 
#                 as.matrix(obs[,c("x","y")]),
#                 as.matrix(cbind(rep(50/nrow(obs),nrow(obs)),obs$cov)),
#                 obs$sett)
# equal proportion
# bs <- lcube(rep(sz/nrow(obs),nrow(obs)), 
#                 as.matrix(obs[,c("x","y")]),
#                 as.matrix(cbind(rep(sz/nrow(obs),nrow(obs)),obs$sett,obs$cov))
#             )