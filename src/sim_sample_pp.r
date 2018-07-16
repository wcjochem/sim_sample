
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
library(splitstackshape)
library(Metrics)
library(parallel)

dir <- "/home/wcj1n15/sim_sample"

## load code to produce model-based estimates ##
source(paste0(dir,"/model_est.r"))
INLA:::inla.dynload.workaround()

# Function to draw samples from a multivariate normal distribution
rmvn <- function(n, mu = 0, V = matrix(1)) {
  p <- length(mu)
  if (any(is.na(match(dim(V), p)))) 
    stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(n * p), ncol = p) %*% D + rep(mu, rep(n, p)))
}

# Function to get a small extent sub-region
make_subregion <- function(r, c, csz, regionsize){ # get the size in number of cells for region
  rcols <- lcols <- trows <- brows <- csz # default square sub-region
  # logic checks to prevent subdomain extending beyond border
  if(c+csz > regionsize){
    rcols <- regionsize - c
    lcols <- lcols + (csz - rcols)
  }
  if(c-csz <= 0){
    lcols <- c - 1
    rcols <- csz + (csz-lcols)
  }
  # repeat checks for rows
  if(r+csz > regionsize){
    brows <- regionsize - r
    trows <- csz + (csz - brows)
  }
  if(r-csz <= 0){
    trows <- r - 1
    brows <- csz + (csz - trows)
  }
  # create sub-region extent
  return(extent(count, r-trows, r+brows, c-lcols, c+rcols))
}

# simulation parameters
# nsims <- 4 # number of spatial fields to repeat
regionsize <- 50 # dim of square region
phi <- c(0.5, 0.1, 0.05, 0.001)
ndraws <- 100 # number of repeated samples
nsims <- length(phi)

# implement sampling strategies
# different strategies -- CHANGE HERE
samps <- c("srs","sys","pwgt","pwgt_ovr") #"strs","areawt"
mods <- c("mbg","brt") #,"rf"
# model + data strategies
strats <- paste(rep(mods, each=length(samps)), samps, sep="_")
strats <- c("pr_srs", strats)
# clean labels - for plotting
cleanlabel <- data.frame(strat=c("pr_srs","mbg_srs","pr_strs","pr_areawt",
                                 "mbg_sys","mbg_pwgt","mbg_pwgt_ovr",
                                 "brt_srs","brt_sys","brt_pwgt","brt_pwgt_ovr",
                                 "rf_srs","rf_sys","rf_pwgt","rf_pwgt_ovr"),
                         name=c("SRS","MBG-SRS","Strat RS","Area wgt Strata",
                                "MBG-SYS","MBG-PPS","MBG-PPS+OS",
                                "BRT-SRS","BRT-SYS","BRT-PPS","BRT-PPS+OS",
                                "RF-SRS","RF-SYS","RF-PPS","RF-PPS+OS"),
                         stringsAsFactors=F)
# different sample sizes
# sampsz <- c(50,100,150,200,250,300,350,400) # CHANGE HERE
sampsz <- c(50,100,150,200) # CHANGE HERE

# abundance model parameters
beta0 <- .4 # intercept
beta1 <- 3 # elevation
beta2 <- -1.5 # quadratic
beta3 <- -1 # trend

## create covariate layers
# elevation
data("volcano") # R dataset for example
# convert matrix to raster
# clipping original down to size
elev <- raster(volcano[10:59,5:54],
               xmn=0, xmx=regionsize, ymn=0, ymx=regionsize)
elev <- scale(elev) # z-scores
  # plot(elev)
  
# trend surface - smoothly varying over space
trend <- raster(elev) # copy elevation raster extent
# get cells to calculate surface values
locns <- xyFromCell(trend, 1:ncell(trend))
# calculate smoothly varying value by x,y coords
vals <- -0.01*locns[,'x'] - 0.01*locns[,'y'] + 0.01*locns[,'x']*locns[,'y']
# update raster values
values(trend) <- vals
trend <- scale(trend) # z-scores
  # plot(trend)

## create spatial fields
# pre-allocate storage
countfields <- vector("list", length=length(phi))
# loop to create abundance datasets
# NOTE this could be slow for large nsims and/or regionsize
for(i in 1:length(phi)){
  print(i)
  set.seed(1)
  # Set up a square study region
  simgrid <- expand.grid(1:regionsize, 1:regionsize)
  n <- nrow(simgrid)
  # Set up distance matrix
  distance <- as.matrix(dist(simgrid))
  # Generate random variable
  X <- rmvn(1, rep(0, n), exp(-phi[i] * distance))
  
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

##
# make mesh for spatial models
sim_mesh <- makemesh(bound=as(extent(elev), "SpatialPolygons"))
  # plot(sim_mesh)


print(Sys.time())
# make cluster
cl <- makeCluster(min(nsims, detectCores()-1))
# load materials on nodes
clusterEvalQ(cl, library(INLA))
clusterEvalQ(cl, library(raster))
clusterEvalQ(cl, library(gstat))
clusterEvalQ(cl, library(Metrics))
clusterEvalQ(cl, library(splitstackshape))

clusterEvalQ(cl, source("/home/wcj1n15/sim_sample/model_est.r"))
clusterEvalQ(cl, INLA:::inla.dynload.workaround())


clusterExport(cl, varlist=c("dir","sim_mesh","countfields","elev","trend","regionsize","samps","strats","sampsz","ndraws"))
clusterExport(cl, varlist=c("make_subregion")) # export functions
# establish node ID
# clusterApply(cl, seq_along(cl), function(x) workerID <<- x)

# Loop and evaluate 
outlist <- clusterApply(cl, 1:nsims, function(i){
  conn <- file( sprintf("/home/wcj1n15/sim_sample/out/process_%d.txt" , Sys.getpid()) , open = "a" )
  r <- 1 # counter to fill in the output
  print(i)
	cat(paste0("sim: ", i), sep="\n", file=conn)
  # storage
  # labels of repetitions, sample size, and observed values
  reslabels <- matrix(data=cbind(rep(sampsz, each=ndraws), NA, NA, NA), 
                      nrow=length(sampsz)*ndraws, ncol=4)
  # storage for per pixel errors
  pred_errs <- matrix(NA, nrow=length(sampsz)*ndraws, ncol=length(strats)*2)
  # storage for total/sub population comparison
  pred_pop <- matrix(NA, nrow=length(sampsz)*ndraws, ncol=length(strats)*2)

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
  # get cell numbers for linking to other sub-regions
  domain$cnum <- cellFromXY(settle, domain[,c("x","y")])
  # keep only settled areas b/c conversion takes all cells
  domain <- domain[domain$counts>0,]
  # extract observed (simulated) population values for all locations
  domain$counts <- extract(count, domain[,c("x","y")])
  # "true" average population per pixel
    mean(domain$counts)
    
  # create 'supergrid' for spatial oversample
  settle_agg <- aggregate(settle, fact=10, sum) # change fact for different scales
  # remove cells with no settlement (not in domain)
  settle_agg[settle_agg==0] <- NA
  # make unique cell ID
  settle_agg[!is.na(settle_agg)] <- 1:length(settle_agg[!is.na(settle_agg)])
  # add this info to the domain
  domain$oversamp <- extract(settle_agg, domain[,c("x","y")])
    
  # extract covariates - for model-based estimates
  domain$elev <- extract(elev, domain[,c("x","y")])
  domain$trend <- extract(trend, domain[,c("x","y")])
      
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
    # plot(count); plot(subreg, add=T); points(domain[which.max(domain$counts),c("x","y")]) # example plot
  # identify sub-region within the domain
  subcells <- cellsFromExtent(count, subreg) # find the cell numbers (matches domain and settle and count)
  domain[domain$cnum %in% subcells, "subdomain"] <- 1 # cell is within the extent of sub-region
  domain[is.na(domain$subdomain), "subdomain"] <- 0 # reclass NA for quick subsetting
  # total "true" population in sub-region
  subpop <- sum(domain[domain$subdomain==1, "counts"])
  
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
  
  ### loop over samples, make repeated draws ###
  for(sz in sampsz){ # vary total sample size
    print(sz)
	cat(paste0(" : ", sz), sep="\n", file=conn)

  # repeated sample realisations
    for(d in 1:ndraws){ 
	cat(paste0("  : ", d), sep="\n", file=conn)

      # store sample size as proportion of settled area
      reslabels[r, 2] <- sz / totsettle
      # store "true" populations
      reslabels[r, 3] <- totpop
      reslabels[r, 4] <- subpop
      # set.seed(d) # ?
    # sampling methods:
    ## simple random sample ##
      if("pr_srs" %in% strats){
        srs <- sample(1:nrow(domain), size=sz, replace=F)
        # extract values of sampled points
        srs <- domain[srs,]
          # plot(count); points(srs[,c("x","y")]) # sample locns
        # sample mean pop per pixel
        samp_mean <- mean(srs$counts)
        # apply mean to settled area domain
        domain$pr_srs <- samp_mean
      }
      
    ## Model-based estimates from random sample ##
    # Geostats #
      if("mbg_srs" %in% strats){
        mbg_srs <- mbg(samp=srs,
                       pred=domain,
                       bound=as(extent(count), "SpatialPolygons"),
                       mesh=sim_mesh)
        # extract the predictions for each location 
        domain$mbg_srs <- mbg_srs$predvals
      }

    ## stratified random sample - equal weight ##
      # equal sample size per stratum
      if("pr_strs" %in% strats){
        szs <- round(sz / nstrat)
        # draw sample
        strs <- data.frame(stratified(domain, "simp_strat", szs))
          # plot(count); points(strs[,c("x","y")])
        # mean values of counts per stratum
        strs_mean <- aggregate(list("mean"=strs$counts), by=list("strat"=strs$simp_strat), FUN=mean)
        # apply mean to settled area domain by stratum
        domain$pr_strs <- strs_mean[domain$simp_strat, "mean"] # expand by domain
      }

   ## stratified random sample - prop to settled area ##
      if("pr_areawt" %in% strats){
        wgt <- zonal(settle, strat, sum)
        wgt <- wgt[,'value'] / totsettle
        # sample sizes - note rounding may causes some different totals
        szs <- setNames(round(wgt*sz), unique(values(strat)))
        # draw sample
        strs <- data.frame(stratified(domain, "simp_strat", szs))
        # mean values of counts per stratum
        strs_mean <- aggregate(list("mean"=strs$counts), by=list("strat"=strs$simp_strat), FUN=mean)
        # overall mean corrected for selection weights
        # sum(wgt * strs_mean$mean)
        # apply mean to settled area domain by stratum
        domain$pr_areawt <- strs_mean[domain$simp_strat, "mean"]        
      }
      
   ## systematic random sample - using oversample spaces as strata ##
      if("mbg_sys" %in% strats){
        # sample size - not rounding may cause some different totals
        szs <- round(sz / length(unique(domain$oversamp)))
        # random sample within systematic strata
        sys <- data.frame(stratified(domain, "oversamp", szs))
        # check sizes and make up for dropped sampled
        # caused by rounding and by some domain areas having too few cells
        miss <- sz - nrow(sys)
        if(miss > 0){#
          avail <- domain[!domain$cnum %in% sys$cnum,] # available locns
          extra_samp <- sample(1:nrow(avail), size=miss, replace=F)
          # add to main sample
          sys <- rbind(sys, avail[extra_samp,])
        }
        # model results
        mbg_sys <- mbg(samp=sys,
                       pred=domain,
                       bound=as(extent(count), "SpatialPolygons"),
                       mesh=sim_mesh)
        # extract predictions
        domain$mbg_sys <- mbg_sys$predvals
      }
      
   ## sample weighted by approximate population density -- test preferential sampling corrections
   ## Model-based estimates ##
      if("mbg_pwgt" %in% strats){
        pwgt_samps <- sample(1:nrow(domain), size=sz, replace=F, prob=domain$pop_wgt)
        # extract values from sampled points
        srs_pwgt <- domain[pwgt_samps,]
        # weighted pop total
        # 1/(srs_pwgt$pop_wgt/sum(domain$pop_wgt)) * (sum(srs_pwgt$counts)/sz)
        # sum(1/(srs_pwgt$pop_wgt/sum(domain$pop_wgt)) * srs_pwgt$counts)/sz
        # (sum(domain$pop_wgt)/srs_pwgt$pop_wgt)/sz * srs_pwgt$counts
        # weighted.mean(srs_pwgt$counts, (sum(domain$pop_wgt)/srs_pwgt$pop_wgt)/sz)
        # model results
        mbg_pwgt <- mbg(samp=srs_pwgt,
                       pred=domain,
                       bound=as(extent(count), "SpatialPolygons"),
                       mesh=sim_mesh)
        # extract the predictions for each location 
        domain$mbg_pwgt <- mbg_pwgt$predvals
      }  
      
    ## spatial oversample using population-weighted sample ##
    ## model-based estimates ##
      if("mbg_pwgt_ovr" %in% strats){
        if(!"mbg_pwgt" %in% strats){ # in case someone skipped previous method
          pwgt_samps <- sample(1:nrow(domain), size=sz, replace=F, prob=domain$pop_wgt)
          # extract values from sampled points
          srs_pwgt <- domain[pwgt_samps,]
        } else{
          srs_pwgt <- domain[pwgt_samps,] # update sample with all domain fields
        }
        # check distribution of sample
        sample_dist <- table(srs_pwgt$oversamp)
        # which areas have no samples
        missed <- !settle_agg[!is.na(settle_agg)] %in% as.numeric(names(sample_dist))
        missed <- settle_agg[!is.na(settle_agg)][missed]
        if(length(missed)>0){
          # find available to drop (must >1 sample in cell)
          avail <- as.numeric(names(sample_dist[sample_dist>1]))
          # drop from sample in available cells a random selection of sites
          drops <- sample(1:nrow(srs_pwgt[srs_pwgt$oversamp %in% avail,]), size=length(missed), replace=F)
          srs_pwgt <- srs_pwgt[-drops,]
          # update with new sample - select 1 per unsampled super cell
          for(c in missed){
            # only select within the large super cell
            newdomain <- domain[domain$oversamp==c,]
            # draw 1 sample in the area and add to the sample
            newsample <- sample(1:nrow(newdomain), size=1, 
                                replace=F, prob=domain[domain$oversamp==c,"pop_wgt"]) # using population weights
            srs_pwgt <- rbind(srs_pwgt, newdomain[newsample,])
          }
	}
        # estimation
        # model results
        mbg_pwgt <- mbg(samp=srs_pwgt,
                       pred=domain,
                       bound=as(extent(count), "SpatialPolygons"),
                       mesh=sim_mesh)
        # extract the predictions for each location 
        domain$mbg_pwgt_ovr <- mbg_pwgt$predvals
      }        
      
   #### Process results ####
      # per-pixel error metrics
      for(s in strats){ # each estimation technique
        idx <- which(strats==s)
        pred_errs[r, idx] <- mape(domain$counts, domain[[s]])
        pred_errs[r, idx+length(strats)] <- rmse(domain$counts, domain[[s]])
      }
      # calculate and store the total and subdomain populations
      for(s in strats){ # each estimation technique
        idx <- which(strats==s)
        pred_pop[r, idx] <- sum(domain[[s]]) # total population
        pred_pop[r, idx+length(strats)] <- sum(domain[domain$subdomain==1, s]) # small area
      }

      r <- r+1
      gc()
    }
  }
  close(conn)
  return(list(i, reslabels, pred_errs, pred_pop))
})
print(Sys.time())
stopCluster(cl)

# save output
saveRDS(outlist, file=paste0(dir,"/out/reslist.rds"))

# End
