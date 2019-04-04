#
# Random sampling from a distributin - 1D space example
#
# Demonstration that with simple random sampling of a non-uniform distribution
# there is a good chance that a small sample will more often come from the "tails."
# The result will be an underestimate of the mean of the distribution.
#
# Chris Jochem, Eric Weber
# 24 May 2018; Update 01 June 2018
#


require(ggplot2)
require(gridExtra)

set.seed(2305)
# repeateed sampling parameters
szs <- c(3,5,10,15,20,25,50) # sample sizes
nsamples <- 1000 # number of repeated samples

# distribution parameters
mu <- 0
sds <- c(0.5, 1, 2, 5, 7)
x <- seq(-10, 10, length=100) # sample locations

# preallocate results storage
reslist <- vector("list", length=length(sds))

# do the loops
for(std in sds){
  # normal distribution to sim 1D region
  y <- dnorm(x, mu, std)
  y_bar <- mean(y) # average value across domain
  
  # storage arrays
  samps <- array(data=NA, dim=c(nsamples, length(szs)))
  sample.means <- array(data=NA, dim=c(nsamples, length(szs))) # array for the sample means
  sample.vars1 <- array(data=NA, dim=c(nsamples, length(szs)))

  # repeat for different size samples
  for(sz in szs){
    for(i in 1:nsamples){
      # random sample
      s <- sample(1:length(x), size=sz, replace=F)
      ys_bar <- mean(y[s]) # sample mean
      
      samps[i, which(szs==sz)] <- mean(y[s] < y_bar) # what proportion are in tails?
      sample.means[i, which(szs==sz)] <- ys_bar # store the sample mean
      sample.vars1[i, which(szs==sz)] <- sum((y[s] - ys_bar)^2) / (sz * (sz-1)) # store variance
    }
  }
  reslist[[which(sds==std)]] <- list(y, y_bar, mean(y < y_bar), 
                                     cbind(1:nsamples,samps),
                                     cbind(1:nsamples,sample.means),
                                     cbind(1:nsamples,sample.vars1))
}

# plots of sample distribution and mean
gplots <- lapply(1:length(reslist), function(j){
  y <- reslist[[j]][[1]]
  y_bar <- reslist[[j]][[2]]
  prop_under <- reslist[[j]][[3]]
  sampdf <- data.frame(reslist[[j]][[4]])
  names(sampdf) <- c("id", paste0("sz_",szs))
  
  sdf.l <- reshape(sampdf,
                   varying=paste0("sz_",szs),
                   v.names="pr_under",
                   timevar="samp_sz",
                   times=szs,
                   direction="long")
  
  p1 <- ggplot(cbind.data.frame(x,y), aes(x=x, y=y)) +
        geom_line(size=1) +
        geom_hline(yintercept=y_bar, col="blue", lty=2, size=1) +
        xlab("1D spatial domain") +
        # geom_text(aes(min(x)+1, y_bar, label="mean", vjust = -1)) +
        theme_bw()
  
  p2 <- ggplot(data=sdf.l, aes(x=as.factor(samp_sz), y=pr_under)) +
        geom_boxplot(fill='#A4A4A4', color="black", show.legend=F) +
        geom_hline(yintercept=prop_under, col="red", size=1) +
        ylab("Proportion low-density sites") +
        xlab("Sample size") +
        theme_bw()
  
  samp.mean.df <- data.frame(reslist[[j]][[5]])
  names(samp.mean.df) <- c("id", paste0("sz_",szs))
  samp.mean.df.reshaped <- reshape(samp.mean.df,
                                   varying=paste0("sz_",szs),
                                   v.names="samp_mean",
                                   timevar="samp_sz",
                                   times=szs,
                                   direction="long")
  
  p3 <- ggplot(data=samp.mean.df.reshaped, aes(x=as.factor(samp_sz), y=samp_mean)) +
        geom_boxplot(fill='gray88', color="gray40", show.legend=F) +
        geom_hline(yintercept=y_bar, col="blue", size=1) +
        stat_summary(fun.y=mean, color="black", fill="green2", geom="point", shape=21, size=3) +
        ylab("Sample means") +
        xlab("Sample size") +
        theme_bw()
  
  vardf <- data.frame(reslist[[j]][[6]])
  names(vardf) <- c("id", paste0("sz_",szs))
  
  vardf.l <- reshape(vardf,
                     varying=paste0("sz_",szs),
                     v.names="samp_var",
                     timevar="samp_sz",
                     times=szs,
                     direction="long")
  # combine mean/variance
  samp_mean_var <- merge(samp.mean.df.reshaped, vardf.l, by=c("id","samp_sz"))
  # upper/lower intervals
  t.95 = qt(0.975, samp_mean_var$samp_sz - 1)
  samp_mean_var$margin <- t.95 * sqrt(samp_mean_var$samp_var)
  samp_mean_var$upper <- samp_mean_var$samp_mean + samp_mean_var$margin
  samp_mean_var$lower <- samp_mean_var$samp_mean - samp_mean_var$margin
  # sort by mean
  samp_mean_var <- samp_mean_var[order(samp_mean_var$samp_sz,samp_mean_var$samp_mean, decreasing=F),]
  samp_mean_var$newid <- rep(1:nsamples, times=length(szs))
  # plot
  p4 <- ggplot(samp_mean_var, aes(x=newid, y=samp_mean, col=as.factor(samp_sz))) +
    geom_linerange(aes(ymin=lower, ymax=upper), col="black", alpha=0.1, show.legend=F) +
    geom_point(show.legend=F) +
    facet_wrap(~ as.factor(samp_sz), scales="free") +
    geom_hline(aes(yintercept=y_bar, group=as.factor(samp_sz)), col="blue", size=1) +
    theme_bw()

  # grid.arrange(p1, p2)
  # plot(arrangeGrob(p1, p2, ncol=2))
  arrangeGrob(p1, p2, p3, ncol=3)
})

# plot the list together
do.call(grid.arrange, c(gplots, ncol=1))




