#
# Process cluster output from sim_sample
# Make graphics
#
# Chris Jochem
# 19 July 2018
#

library(ggplot2)
library(rasterVis)
library(gridExtra)
library(grid)

####
# Function to draw samples from a multivariate normal distribution
rmvn <- function(n, mu = 0, V = matrix(1)) {
  p <- length(mu)
  if (any(is.na(match(dim(V), p)))) 
    stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(n * p), ncol = p) %*% D + rep(mu, rep(n, p)))
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
# make plots of simulated fields
count_plots <- lapply(1:length(phi), function(x){
  f <- countfields[[x]]
  f[f==0] <- NA
  gplot(f) +
    geom_tile(aes(fill= value)) + 
    # scale_fill_gradient(low="white", high="blue") +
    scale_fill_gradientn(colours = terrain.colors(10), na.value="transparent") +
    theme_bw() +
    coord_equal()
})
# arrange all plots
marrangeGrob(count_plots, nrow=2, ncol=2)



#####
# process cluster output
reslist <- readRDS("C:/Users/wcj1n15.SOTON/Dropbox/Soton/proj/Nigeria/density_all/sim/reslist_2018-07-24.rds")
  length(reslist)
# list of lists: i, reslabels, pred_errs, pred_pop
errs <- vector("list", length=length(reslist))
preds <- vector("list", length=length(reslist))
#
for(i in 1:length(reslist)){
  l <- reslist[[i]]
  n <- l[[1]]
  labs <- cbind(n, l[[2]])
  errs[[i]] <- cbind(labs, l[[3]])
  preds[[i]] <- cbind(labs, l[[4]])
}
# combine lists
pred_errs <- do.call(rbind, errs)
pred_pop <- do.call(rbind, preds)


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
  facet_wrap(~sim, scales="free", ncol=4) +
  ggtitle("RMSE") +
  ylab("RMSE") +
  xlab("Sample size") +
  theme_bw()

ylim1 = boxplot.stats(err.rmse.l$err)$stats[c(1, 5)]
# scale y limits based on ylim1
grmse <- grmse + coord_cartesian(ylim = ylim1*1.1)

gmape <- ggplot(data=err.mape.l, aes(x=as.factor(sz), y=err, fill=est)) +
  geom_boxplot() +
  facet_wrap(~sim, scales="free", ncol=4) +
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
  facet_wrap(~sim, scales="free", nrow=4) +
  geom_hline(data=totpop, aes(yintercept=totpop, col="red"), show.legend=F) +
  ggtitle("Estimated total population") +
  ylab("Population") +
  xlab("Sample size") +
  theme_bw() +
  theme(legend.position="bottom")

# sub regions
hisubpop.df.l <- hisubpop.df.l[!is.na(hisubpop.df.l$pop) & !is.infinite(hisubpop.df.l$pop),]
hisubpop.df.l <- subset(hisubpop.df.l, pop < hisubpop*3)

losubpop.df.l <- losubpop.df.l[!is.na(losubpop.df.l$pop) & !is.infinite(losubpop.df.l$pop),]
losubpop.df.l <- subset(losubpop.df.l, pop < losubpop*3)

# plot prediction of subregion population
ghisubpop <- ggplot(data=hisubpop.df.l, aes(x=as.factor(sz), y=pop, fill=est)) + 
  geom_boxplot() +
  facet_wrap(~sim, scales="free", nrow=4) +
  geom_hline(data=totsubpop_hi, aes(yintercept=hisubpop, col="red"), show.legend=F) +
  ggtitle("Estimated high population subregion") +
  ylab("Population") +
  xlab("Sample size") +
  theme_bw() +
  theme(legend.position="none")

glosubpop <- ggplot(data=losubpop.df.l, aes(x=as.factor(sz), y=pop, fill=est)) + 
  geom_boxplot() +
  facet_wrap(~sim, scales="free", nrow=4) +
  geom_hline(data=totsubpop_lo, aes(yintercept=losubpop, col="red"), show.legend=F) +
  ggtitle("Estimated low population subregion") +
  ylab("Population") +
  xlab("Sample size") +
  theme_bw() +
  theme(legend.position="none")

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
  # return(combined)
}  

# grid_arrange_shared_legend(gtotpop, gsubpop, grmse, gmape, ncol=4, nrow=1)
grid_arrange_shared_legend(gtotpop, ghisubpop, glosubpop, ncol=3, nrow=1)
grid_arrange_shared_legend(grmse, gmape, ncol=1, nrow=2)

##
# library(cowplot)
# 
# # get legend of boxplots
# legend <- get_legend(gtotpop)
# gtotpopa <- gtotpop + theme(legend.position = "none")
# # plot_grid(gtotpopa, ghisubpop, glosubpop, nrow=1, ncol=3)
# 
# # component plots
# # simulations
# simplots <- align_plots(plotlist=count_plots, align='v', axis='l')
# allsim <- plot_grid(plotlist=simplots, nrow=4)
# 
# # align all boxplots vertically
# popplots <- align_plots(gtotpopa, ghisubpop, glosubpop, align='v', axis='l')
# toprow <- plot_grid(popplots[[1]],popplots[[2]],popplots[[3]], ncol=3)
# # combined with legend
# pop <- plot_grid(toprow, legend, nrow=2, rel_heights=c(1,.1))
# 
# # put together everything
# plot_grid(allsim, pop, ncol=2)
# 
# simpop <- plot_grid(allsim, toprow, ncol=2)
# 
# 
