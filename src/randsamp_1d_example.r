#
# Random sampling from a distributin - 1D space example
#
# Demonstration that with simple random sampling of a non-uniform distribution
# there is a good chance that a small sample will more often come from the "tails."
# The result will be an underestimate of the mean of the distribution.
#
# Chris Jochem
# 24 May 2018
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
  
  # storage array
  samps <- array(data=NA, dim=c(nsamples, length(szs)))
  # repeat for different size samples
  for(sz in szs){
    for(i in 1:nsamples){
      # random sample
      s <- sample(1:length(x), size=sz, replace=F)
      ys_bar <- mean(y[s]) # sample mean
      
      samps[i, which(szs==sz)] <- mean(y[s] < y_bar) # what proportion are in tails?
    }
  }
  reslist[[which(sds==std)]] <- list(y, y_bar, mean(y < y_bar), cbind(1:nsamples,samps))
}

# plots
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
  
  # grid.arrange(p1, p2)
  # plot(arrangeGrob(p1, p2, ncol=2))
  arrangeGrob(p1, p2, ncol=2)
})

# plot the list together
do.call(grid.arrange, c(gplots, ncol=1))
