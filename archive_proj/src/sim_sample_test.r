
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
library(grid)
library(splitstackshape)
library(Metrics)

## load code to produce model-based estimates ##
source("./src/model_est.r")

## HELPER FUNCTIONS ##
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
# phi <- 0.05 # smoothness parameter
phi <- c(0.5, 0.1, 0.05, 0.001)
ndraws <- 10 # number of repeated samples
nsims <- length(phi)

# implement sampling strategies
# different strategies -- CHANGE HERE
# list used to automate data storage creation and logic to limit evaluation steps
samps <- c("srs","sys","pwgt","pwgt_ovr") #"strs","areawt"
mods <- c("mbg","brt","rf") #,"rf"
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
sampsz <- c(50)

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
countfields <- vector("list", length=length(phi))
# loop to create abundance datasets
# NOTE this could be slow for large nsims and/or regionsize
for(i in 1:length(phi)){
  print(i)
  set.seed(1) # fixed seed - varying phi
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

##
  par(mfrow=c(2,2), mar=c(1,1,1,1))
  plot(countfields[[1]], main="Sim 1")
  plot(countfields[[2]], main="Sim 2")
  plot(countfields[[3]], main="Sim 3")
  plot(countfields[[4]], main="Sim 4")
  par(mfrow=c(1,1))

## add additional simulated surfaces ##
# # homongenous 
# counts <- 4
# homog_pop <-rasterFromXYZ(cbind(coordinates(elev), counts))

# # purely random
# counts <- rnorm(ncell(elev), mean=5, sd=.5) # random values
# rand_pop <- rasterFromXYZ(cbind(coordinates(elev), counts))
# 
# # purely covariate
# lambda <- exp(beta0 + beta1*values(elev) + beta2*values(elev)^2 + beta3*values(trend)) # same generating function
# # generate counts from the simulated intensity
# counts <- rpois(n, lambda)
#     # sum(counts)
# # gridded abundance counts
# cov_pop <- rasterFromXYZ(cbind(coordinates(elev), counts))
#   # plot(cov_pop)
# # combine with list of spatial field simulations
# addfields <- list(rand_pop, cov_pop)
# # countfields <- c(addfields, countfields)
# # update the number of simulations
# nsims <- length(countfields)

# visualise the results
#   plot(countfields[[1]]) # first example
# create list of all counts
count_plots <- lapply(1:length(phi), function(x){
  gplot(countfields[[x]]) +
    geom_tile(aes(fill= value)) + 
    scale_fill_gradient(low="white", high="red") + 
    theme_bw() +
    coord_equal()
})
  # count_plots
# arrange all plots
marrangeGrob(count_plots, nrow=2, ncol=2)


##
# preallocate storage for the output
# labels of repetitions, sample size, and observed values
reslabels <- matrix(data=cbind(rep(1:nsims, each=length(sampsz)*ndraws), 
                               rep(sampsz, each=ndraws), NA, NA, NA, NA), 
                    nrow=nsims*length(sampsz)*ndraws, ncol=6)
# storage for per pixel errors
pred_errs <- matrix(NA, nrow=nsims*length(sampsz)*ndraws, ncol=length(strats)*2) # x2 for RMSE and MAPE
# storage for total/sub population comparison
pred_pop <- matrix(NA, nrow=nsims*length(sampsz)*ndraws, ncol=length(strats)*3) # x3 for total, hi/low domains
# storage to record sample locations from SRS
srs_locn <- vector("list", length=nsims)

# make mesh for spatial models
sim_mesh <- makemesh(bound=as(extent(elev), "SpatialPolygons"))
  # plot(sim_mesh)

# Loop and evaluate 
r <- 1 # counter to fill in the output
print(Sys.time())
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
    
  # store the domain for sample counts
  srs_locn[[i]] <- domain
    
  # create a high pop sub-domain area to estimate totals
  # analogous to having a capital city within a regional sampling domain
  hicell <- cellFromXY(count, domain[which.max(domain$counts),c("x","y")])
  row <- rowFromCell(count, hicell) # get row/column of the high population area
  col <- colFromCell(count, hicell)

  csz <- 4 # set sub-area size (square: csz+1 x csz+1 cells)
  stopifnot(regionsize >= 2*csz + 1) # double-check no small regions
  # create sub-region
  subreg <- make_subregion(row, col, csz, regionsize)
    plot(count); plot(subreg, add=T); points(domain[which.max(domain$counts),c("x","y")]) # example plot
  # identify sub-region within the domain
  subcells <- cellsFromExtent(count, subreg) # find the cell numbers (matches domain and settle and count)
  domain[domain$cnum %in% subcells, "hidomain"] <- 1 # cell is within the extent of sub-region
  domain[is.na(domain$hidomain), "hidomain"] <- 0 # reclass NA for quick subsetting
  # total "true" population in sub-region
  hisubpop <- sum(domain[domain$hidomain==1, "counts"])
  
  # create a low pop sub-domain area to estimate totals
  locell <- cellFromXY(count, domain[which.min(domain$counts),c("x","y")])
  row <- rowFromCell(count, locell) # get row/column of the high population area
  col <- colFromCell(count, locell)
  
  csz <- 4 # set sub-area size (square: csz+1 x csz+1 cells)
  stopifnot(regionsize >= 2*csz + 1) # double-check no small regions
  # create sub-region
  subreg <- make_subregion(row, col, csz, regionsize)
    # plot(count); plot(subreg, add=T); points(domain[which.min(domain$counts),c("x","y")]) # example plot
  # identify sub-region within the domain
  subcells <- cellsFromExtent(count, subreg) # find the cell numbers (matches domain and settle and count)
  domain[domain$cnum %in% subcells, "lodomain"] <- 1 # cell is within the extent of sub-region
  domain[is.na(domain$lodomain), "lodomain"] <- 0 # reclass NA for quick subsetting
  # total "true" population in sub-region
  losubpop <- sum(domain[domain$lodomain==1, "counts"])
  
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
    # add column for each sample size
    srs_locn[[i]][[paste0("sz_",sz)]] <- 0

  # repeated sample realisations
    for(d in 1:ndraws){
      print(paste0(" ", d))
      # store sample size as proportion of settled area
      reslabels[r, 3] <- sz / totsettle
      # store "true" populations
      reslabels[r, 4] <- totpop
      reslabels[r, 5] <- hisubpop
      reslabels[r, 6] <- losubpop
      # set.seed(d) # ?
    # sampling methods:
    ## simple random sample ##
      if("srs" %in% samps){
        srs <- sample(1:nrow(domain), size=sz, replace=F)
        # extract values of sampled points
        srs <- domain[srs,]
          # plot(count); points(srs[,c("x","y")]) # sample locns
        # store iterative count of selected locns
        srs_locn[[i]][srs_locn[[i]]$cnum %in% srs$cnum, paste0("sz_",sz)] <- 
          srs_locn[[i]][srs_locn[[i]]$cnum %in% srs$cnum, paste0("sz_",sz)] + 1
        
        # simple mean
        if("pr_srs" %in% strats){
          # sample mean pop per pixel
          samp_mean <- mean(srs$counts)
          # apply mean to settled area domain
          domain$pr_srs <- samp_mean
        }
        # geostats
        if("mbg_srs" %in% strats){
          mbg_srs <- mbg(samp=srs,
                       pred=domain,
                       bound=as(extent(count), "SpatialPolygons"),
                       mesh=sim_mesh)
          # extract the predictions for each location 
          domain$mbg_srs <- mbg_srs$predvals
        }
        # boosted regression tree
        if("brt_srs" %in% strats){
          brt_srs <- brt(samp=srs,
                         pred=domain)
          # extraction predictions
          domain$brt_srs <- brt_srs$predvals
        }
        # random forest
        if("rf_srs" %in% strats){
          rf_srs <- rf(samp=srs,
                       pred=domain)
          # extract predictions
          domain$rf_srs <- rf_srs$predvals
        }
      }

    ## systematic random sample - using oversample spaces as strata ##
      if("sys" %in% samps){
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
        # geostats
        if("mbg_sys" %in% strats){
          mbg_sys <- mbg(samp=sys,
                         pred=domain,
                         bound=as(extent(count), "SpatialPolygons"),
                         mesh=sim_mesh)
          # extract predictions
          domain$mbg_sys <- mbg_sys$predvals
        }
        # boosted regression tree
        if("brt_sys" %in% strats){
          brt_sys <- brt(samp=sys,
                         pred=domain)
          # extract predictions
          domain$brt_sys <- brt_sys$predvals
        }
        # random forest
        if("rf_sys" %in% strats){
          rf_sys <- rf(samp=sys,
                       pred=domain)
          # extract predictions
          domain$rf_sys <- rf_sys$predvals
        }
      }  

    ## sample weighted by approximate population density -- test preferential sampling corrections
      if("pwgt" %in% samps){
        pwgt_samps <- sample(1:nrow(domain), size=sz, replace=F, prob=domain$pop_wgt)
        # extract values from sampled points
        srs_pwgt <- domain[pwgt_samps,]
        # weighted pop total
        # 1/(srs_pwgt$pop_wgt/sum(domain$pop_wgt)) * (sum(srs_pwgt$counts)/sz)
        # sum(1/(srs_pwgt$pop_wgt/sum(domain$pop_wgt)) * srs_pwgt$counts)/sz
        # (sum(domain$pop_wgt)/srs_pwgt$pop_wgt)/sz * srs_pwgt$counts
        # weighted.mean(srs_pwgt$counts, (sum(domain$pop_wgt)/srs_pwgt$pop_wgt)/sz)
      
        # geostats
        if("mbg_pwgt" %in% strats){
          mbg_pwgt <- mbg(samp=srs_pwgt,
                          pred=domain,
                          bound=as(extent(count), "SpatialPolygons"),
                          mesh=sim_mesh)
          # extract the predictions for each location 
          domain$mbg_pwgt <- mbg_pwgt$predvals
        }
        # boosted regression tree
        if("brt_pwgt" %in% strats){
          brt_pwgt <- brt(samp=srs_pwgt,
                          pred=domain)
          # extract predictions
          domain$brt_pwgt <- brt_pwgt$predvals
        }
        # random forest
        if("rf_pwgt" %in% strats){
          rf_pwgt <- rf(samp=srs_pwgt,
                        pred=domain)
          # extract predictions
          domain$rf_pwgt <- rf_pwgt$predvals
        }
      }  

    ## spatial oversample using population-weighted sample ##
      if("pwgt_ovr" %in% samps){
        if(!"pwgt" %in% samps){ # in case skipped previous method
          pwgt_samps <- sample(1:nrow(domain), size=sz, replace=F, prob=domain$pop_wgt)
        }
        # extract values from sampled points
        srs_pwgt <- domain[pwgt_samps,] # update sample with all domain fields
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
        # geostats
        if("mbg_pwgt_ovr" %in% strats){
          mbg_pwgt_ovr <- mbg(samp=srs_pwgt,
                              pred=domain,
                              bound=as(extent(count), "SpatialPolygons"),
                              mesh=sim_mesh)
          # extract the predictions for each location 
          domain$mbg_pwgt_ovr <- mbg_pwgt_ovr$predvals
        }
        # boosted regression tree
        if("brt_pwgt_ovr" %in% strats){
          brt_pwgt_ovr <- brt(samp=srs_pwgt,
                              pred=domain)
          # extract predictions
          domain$brt_pwgt_ovr <- brt_pwgt_ovr$predvals
        }
        # random forest
        if("rf_pwgt_ovr" %in% strats){
          rf_pwgt_ovr <- rf(samp=srs_pwgt,
                            pred=domain)
          # extract predictions
          domain$rf_pwgt_ovr <- rf_pwgt_ovr$predvals
        } 
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
        pred_pop[r, idx+length(strats)] <- sum(domain[domain$hidomain==1, s]) # small high area
        pred_pop[r, idx+(2*length(strats))] <- sum(domain[domain$lodomain==1, s]) # small low area
      }

      r <- r+1
      gc()
    }
  }
}
print(Sys.time())

# add the labels to matrix of results
pred_errs <- cbind(reslabels, pred_errs)
pred_pop <- cbind(reslabels, pred_pop)

####
# # process cluster output
# reslist <- readRDS("C:/Users/wcj1n15.SOTON/Dropbox/Soton/proj/Nigeria/density_all/sim/reslist.rds")
#   length(reslist)
# # list of lists: i, reslabels, pred_errs, pred_pop
# errs <- vector("list", length=length(reslist))
# preds <- vector("list", length=length(reslist))
# #
# for(i in 1:length(reslist)){
#   l <- reslist[[i]]
#   n <- l[[1]]
#   labs <- cbind(n, l[[2]])
#   errs[[i]] <- cbind(labs, l[[3]])
#   preds[[i]] <- cbind(labs, l[[4]])
# }
# 
# pred_errs <- do.call(rbind, errs)
# pred_pop <- do.call(rbind, preds)

####

#####
# process results
# convert errors to datafraem for ggplot
err.df <- data.frame(pred_errs)
names(err.df) <- c("sim","sz","pct_samp","totpop","hisubpop","losubpop", paste0("mape_", strats), paste0("rmse_", strats))
# split
err.df.mape <- err.df[,c("sim","sz","pct_samp", paste0("mape_", strats))]
err.df.rmse <- err.df[,c("sim","sz","pct_samp",paste0("rmse_", strats))]
# convert to long format
err.rmse.l <- reshape(err.df.rmse,
                    varying=paste0("rmse_", strats),
                    v.names="err",
                    timevar="est",
                    times=paste0("rmse_", strats),
                    direction="long")

err.mape.l <- reshape(err.df.mape,
                    varying=paste0("mape_", strats),
                    v.names="err",
                    timevar="est",
                    times=paste0("mape_", strats),
                    direction="long")

err.rmse.l$Method <- cleanlabel[match(gsub("rmse_","", err.rmse.l$est, fixed=T), cleanlabel$strat),"name"]

# plot EA-level error metrics
grmse <- ggplot(data=err.rmse.l, aes(x=as.factor(sz), y=err, fill=Method)) +
  geom_boxplot() +
  facet_wrap(~sim, scales="free", nrow=3) +
  ggtitle("RMSE") +
  ylab("RMSE") +
  xlab("Sample size") +
  theme_bw()

ylim1 = boxplot.stats(err.rmse.l$err)$stats[c(1, 5)]
# scale y limits based on ylim1
grmse <- grmse + coord_cartesian(ylim = ylim1*1.1)

gmape <- ggplot(data=err.mape.l, aes(x=as.factor(sz), y=err, fill=est)) +
  geom_boxplot() +
  facet_wrap(~sim, scales="free", nrow=3) +
  ggtitle("MAPE") +
  ylab("MAPE") +
  xlab("Sample size") +
  theme_bw()

ylim1 = boxplot.stats(err.mape.l$err)$stats[c(1, 5)]
# scale y limits based on ylim1
gmape <- gmape + coord_cartesian(ylim = ylim1*1.1)



# convert population to dataframe for ggplot
pop.df <- data.frame(pred_pop)
names(pop.df) <- c("sim","sz","pct_samp","totpop","hisubpop","losubpop", strats, paste0(strats, "_hi"), paste0(strats, "_lo"))
# subdomain data frame
hisubpop.df <- pop.df[,c("sim","sz","pct_samp","hisubpop", paste0(strats, "_hi"))]
losubpop.df <- pop.df[,c("sim","sz","pct_samp","losubpop", paste0(strats, "_lo"))]

pop.df <- pop.df[,!names(pop.df) %in% c(paste0(strats, "_hi"), paste0(strats, "_lo"))]

# pop records convert to long format
pop.df.l <- reshape(pop.df,
                    varying=strats,
                    v.names="pop",
                    timevar="est",
                    times=strats,
                    direction="long")
# add clean labels
pop.df.l$Method <- cleanlabel[match(pop.df.l$est, cleanlabel$strat),"name"]
# long-format for subpop data frames - for plotting
hisubpop.df.l <- reshape(hisubpop.df,
                         varying=paste0(strats, "_hi"),
                         v.names="pop",
                         timevar="est",
                         times=strats,
                         direction="long")

losubpop.df.l <- reshape(losubpop.df,
                         varying=paste0(strats, "_lo"),
                         v.names="pop",
                         timevar="est",
                         times=strats,
                         direction="long")

# find simulation-specific total population
totpop <- unique(pop.df[,c("sim","totpop")])
totsubpop_hi <- unique(hisubpop.df[,c("sim","hisubpop")])
totsubpop_lo <- unique(losubpop.df[,c("sim","losubpop")])

pop.df.l <- pop.df.l[!is.na(pop.df.l$pop) & !is.infinite(pop.df.l$pop),]
pop.df.l <- subset(pop.df.l, pop < totpop*3)

# plot total popuation predictions
gtotpop <- ggplot(data=pop.df.l, aes(x=as.factor(sz), y=pop, fill=Method)) + 
  geom_boxplot() +
  facet_wrap(~sim, scales="free", nrow=3) +
  geom_hline(data=totpop, aes(yintercept=totpop, col="red"), show.legend=F) +
  ggtitle("Estimated total population") +
  ylab("Population") +
  xlab("Sample size") +
  theme_bw()

# sub regions
hisubpop.df.l <- hisubpop.df.l[!is.na(hisubpop.df.l$pop) & !is.infinite(hisubpop.df.l$pop),]
hisubpop.df.l <- subset(hisubpop.df.l, pop < hisubpop*3)

losubpop.df.l <- losubpop.df.l[!is.na(losubpop.df.l$pop) & !is.infinite(losubpop.df.l$pop),]
hisubpop.df.l <- subset(losubpop.df.l, pop < losubpop*3)

# plot prediction of subregion population
ghisubpop <- ggplot(data=hisubpop.df.l, aes(x=as.factor(sz), y=pop, fill=est)) + 
  geom_boxplot() +
  facet_wrap(~sim, scales="free", nrow=3) +
  geom_hline(data=totsubpop_hi, aes(yintercept=hisubpop, col="red"), show.legend=F) +
  ggtitle("Estimated subregion population") +
  ylab("Population") +
  xlab("Sample size") +
  theme_bw()

glosubpop <- ggplot(data=losubpop.df.l, aes(x=as.factor(sz), y=pop, fill=est)) + 
  geom_boxplot() +
  facet_wrap(~sim, scales="free", nrow=3) +
  geom_hline(data=totsubpop_lo, aes(yintercept=losubpop, col="red"), show.legend=F) +
  ggtitle("Estimated subregion population") +
  ylab("Population") +
  xlab("Sample size") +
  theme_bw()

## plot combined results
# res <- arrangeGrob(gtotpop, gsubpop, grmse, gmape, ncol=2)
  # plot(res)
  
####
# https://github.com/tidyverse/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs  
grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)

  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  grid.newpage()
  grid.draw(combined)
  # return gtable invisibly
  invisible(combined)
}  

# grid_arrange_shared_legend(gtotpop, gsubpop, grmse, gmape, ncol=4, nrow=1)
grid_arrange_shared_legend(gtotpop, ghisubpop, glosubpop, ncol=3, nrow=1)
grid_arrange_shared_legend(grmse, gmape, ncol=2, nrow=1)

# End
saveRDS(countfields, file="C:/Users/wcj1n15.SOTON/Dropbox/Soton/proj/mcs_design/countfields.rds")
saveRDS(pred_errs, file="C:/Users/wcj1n15.SOTON/Dropbox/Soton/proj/mcs_design/pred_errs.rds")
saveRDS(pred_pop, file="C:/Users/wcj1n15.SOTON/Dropbox/Soton/proj/mcs_design/pred_pop.rds")
