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

# stratified random by type
szs <- trunc(sz / length(unique(obs$sett)))
sts <- data.frame(stratified(obs, "sett", szs))$rid
# redistribute other samples randomly
rem <- sz %% length(unique(obs$sett))
extra_sts <- sample(obs[!obs$rid %in% sts,"rid"], rem, replace=F)
sts <- c(sts, extra_sts) # combine

# systematic random by admin2
szs <- trunc(sz / length(unique(obs$adm2)))
sys <- data.frame(stratified(obs, "adm2", szs))$rid
# redistribute other samples randomly
rem <- sz %% length(unique(obs$adm2))
# randomly select blocks
extra_adm <- sample(unique(obs$adm2), size=rem, replace=F)
# sample a pixel within admin unit
extra <- data.frame(stratified(obs[!obs$rid %in% sys,], "adm2", 1, select=list(adm2=extra_adm)))$rid
sys <- c(sys,extra) # combine

# spatially-balanced
# population-weighted
bs <- lcube(prob=wts, 
            Xspread=as.matrix(obs[,c("x","y")]),
            Xbal=as.matrix(cbind(wts,obs[,c("cov","sett","x","y")]))
           )

# set samples
obs[srs, "srs"] <- TRUE
obs[is.na(obs$srs),"srs"] <- FALSE
obs[pwgt, "pwt"] <- TRUE
obs[is.na(obs$pwt), "pwt"] <- FALSE
obs[obs$rid %in% sts, "sts"] <- TRUE
obs[is.na(obs$sts), "sts"] <- FALSE
obs[obs$rid %in% sys, "sys"] <- TRUE
obs[is.na(obs$sys), "sys"] <- FALSE
obs[bs, "bal"] <- TRUE
obs[is.na(obs$bal), "bal"] <- FALSE


#####
# Model Fitting
#####
c.c <- list(cpo=TRUE, dic=TRUE, waic=TRUE, config=TRUE)
nsamp <- 1e3 # number of samples for predictions
# create spatial mesh for geostatistical models
fake.locations <- matrix(c(0,0,r_width,r_height,0,r_height,r_width,0), nrow=4, byrow=T)
sim_mesh <- inla.mesh.2d(loc=fake.locations, max.edge=c(1, 4)) 
  # plot(sim_mesh, main="")
# spde - spatial prior
spde <- inla.spde2.pcmatern(sim_mesh, prior.range=c(10, .9), prior.sigma=c(.5, .5))



### Simple Random Sample ###

dat <- obs[obs$srs==T,]
pred_a <- obs[obs$srs==F,]
## Geostats model - Fixed Intercept
# set up model
A.est <- inla.spde.make.A(mesh=sim_mesh,
                          loc=data.matrix(dat[,c("x","y")]))

A.pred_in <- inla.spde.make.A(mesh=sim_mesh,
                              loc=data.matrix(pred_a[,c("x","y")]))

A.pred_out <- inla.spde.make.A(mesh=sim_mesh,
                               loc=data.matrix(pred[,c("x","y")]))
# index to the mesh
mesh.index0 <- inla.spde.make.index(name="field", n.spde=spde$n.spde)


# data stack for model
stack.est <- inla.stack(data=list(pop=dat$pop),
                        A=list(A.est, 1),
                        effects=list( c(mesh.index0, list(Intercept=1)),
                                      list(dat[,c("cov","sett")]) ),
                        tag='est')
# fit model
form <- pop ~ -1 + Intercept + cov + f(sett, model="iid", values=1:3) + f(field, model=spde)

res <- inla(form, 
           family="poisson",
           data=inla.stack.data(stack.est),
           control.predictor=list(A=inla.stack.A(stack.est), compute=T, link=1),
           control.compute=c.c) 
# prediction
# draw samples from the posterior
ps <- inla.posterior.sample(n=nsamp, res)
# get indices to the effects
contents <- res$misc$configs$contents

idSpace <- contents$start[which(contents$tag=="field")]-1 + (1:contents$length[which(contents$tag=="field")])
idSett <- contents$start[which(contents$tag=="sett")]-1 + (1:contents$length[which(contents$tag=="sett")])
idX <- contents$start[which(contents$tag=="Intercept")]-1 + (1:2) # fixed effects, change for covariates
# extract samples
xLatent <- matrix(0, nrow=length(ps[[1]]$latent), ncol=nsamp)
xHyper <- matrix(0, nrow=length(ps[[1]]$hyperpar), ncol=nsamp)
for(i in 1:nsamp){
  xLatent[,i] <- ps[[i]]$latent
  xHyper[,i] <- ps[[i]]$hyperpar
}
xSpace <- xLatent[idSpace,]
xSett <- xLatent[idSett,]
xX <- xLatent[idX,]

# construct predictions
# in-sample
linpred <- as.matrix(A.est %*% xSpace + xSett[dat$sett,] + as.matrix(cbind(1, dat[,c("cov")])) %*% xX)
allpred_a <- cbind(dat[dat$srs==T,], linpred)
pred_N_s <- data.frame(t(apply(linpred, 1, FUN=function(x){ quantile(rpois(n=nsamp, lambda=exp(x)), probs=c(0.025,0.5,0.975)) })))
names(pred_N_s) <- c("lower","median","upper")

# within sample area
linpred <- as.matrix(A.pred_in %*% xSpace + xSett[pred_a$sett,] + as.matrix(cbind(1, pred_a[,c("cov")])) %*% xX)
allpred_a <- rbind(allpred_a, cbind(pred_a, linpred))
pred_N_in <- data.frame(t(apply(linpred, 1, FUN=function(x){ quantile(rpois(n=nsamp, lambda=exp(x)), probs=c(0.025,0.5,0.975)) })))
names(pred_N_in) <- c("lower","median","upper")
# outside sample area
linpred <- as.matrix(A.pred_out %*% xSpace + xSett[pred$sett,] + as.matrix(cbind(1, pred[,c("cov")])) %*% xX)
  quantile(apply(linpred, 2, FUN=sum), probs=c(0.025, 0.5, 0.975))
pred_N_out <- data.frame(t(apply(linpred, 1, FUN=function(x){ quantile(rpois(n=nsamp, lambda=exp(x)), probs=c(0.025,0.5,0.975)) })))
names(pred_N_out) <- c("lower","median","upper")

# plot: obs vs. pred in sample
plotdat <- pred_N_s
plotdat$obs <- dat$pop
plotdat$pred <- plotdat$median

plot(NA,
     main=paste('Population Density \n(r2 =',round(cor(plotdat$pred, plotdat$obs)^2,2),')'),
     ylim=c(0,max(plotdat$upper)),
     xlim=c(0, max(plotdat$obs)),
     xlab='observed', ylab='predicted')
for(i in 1:nrow(plotdat)){
  arrows(x0=plotdat[i,'obs'], y0=plotdat[i,'lower'], y1=plotdat[i,'upper'], length=0, col=rgb(0.5,0.5,0.5,0.5))
}
points(plotdat$obs, plotdat$pred)
abline(0,1,col='red')

# plot: obs vs. pred in domain
plotdat <- pred_N_in
plotdat$obs <- pred_a$pop
plotdat$pred <- plotdat$median

plot(NA,
     main=paste('Population Density \n(r2 =',round(cor(plotdat$pred, plotdat$obs)^2,2),')'),
     ylim=c(0,max(plotdat$upper)),
     xlim=c(0, max(plotdat$obs)),
     xlab='observed', ylab='predicted')
for(i in 1:nrow(plotdat)){
  arrows(x0=plotdat[i,'obs'], y0=plotdat[i,'lower'], y1=plotdat[i,'upper'], length=0, col=rgb(0.5,0.5,0.5,0.5))
}
points(plotdat$obs, plotdat$pred)
abline(0,1,col='red')

# plot: obs vs. pred out of domain
plotdat <- pred_N_out
plotdat$obs <- pred$pop
plotdat$pred <- plotdat$median

plot(NA,
     main=paste('Population Density \n(r2 =',round(cor(plotdat$pred, plotdat$obs)^2,2),')'),
     ylim=c(0,max(plotdat$upper)),
     xlim=c(0, max(plotdat$obs)),
     xlab='observed', ylab='predicted')
for(i in 1:nrow(plotdat)){
  arrows(x0=plotdat[i,'obs'], y0=plotdat[i,'lower'], y1=plotdat[i,'upper'], length=0, col=rgb(0.5,0.5,0.5,0.5))
}
points(plotdat$obs, plotdat$pred)
abline(0,1,col='red')

## "town" populations
quantile(apply(allpred_a[allpred_a$town==1,16:1015], 2, function(x) sum(x, na.rm=T)), probs=c(0.025, 0.5, 0.975))
quantile(apply(allpred_a[allpred_a$town==4,16:1015], 2, function(x) sum(x, na.rm=T)), probs=c(0.025, 0.5, 0.975))

## Non-spatial - Random Intercept



## Geostats model - Random Intercept










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