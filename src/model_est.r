#
# Code for simplified model-based estimates
# Used with sim_sample_test.r 
#
# Chris Jochem
# 05 June 2018
#

require(inlabru) 
require(INLA)

## simplified Bayesian geostatistcal model ##
# make a mesh for SPDE model
makemesh <- function(locn=NULL, bound=NULL, me=c(1,2), co=1, off=c(5,10)){
  mesh <- inla.mesh.2d(loc=locn,
                       boundary=bound,
                       max.edge=me,
                       cutoff=co,
                       offset=off)

  return(mesh)
}


# general function for models
mbg <- function(samp=NULL, pred=NULL, bound=NULL, mesh=NULL,
                meshconfig=list(me=c(1,2), co=1, off=c(5,10)) ){
  
  if(is.null(mesh)){
    mesh <- makemesh(samp[,c("x","y")], bound)
    # plot(mesh); points(samp[,c("x","y")], pch=16, col="red")
  }
  # spde model - penalised prior
  matern <- inla.spde2.pcmatern(mesh,
                                prior.sigma=c(10, 0.01),
                                prior.range=c(15, 0.50))
  
  # components
  form <- l_counts ~ -1 + Intercept + elev + trend + f(field, model=matern)
  c.c <- list(cpo=TRUE, dic=TRUE, waic=TRUE, config=TRUE)
  # set up model
  A.est <- inla.spde.make.A(mesh=mesh,
                            loc=data.matrix(samp[,c("x","y")]))
  
  A.pred <- inla.spde.make.A(mesh=mesh,
                             loc=data.matrix(pred[,c("x","y")]))
  # index to the mesh
  mesh.index0 <- inla.spde.make.index(name="field", n.spde=matern$n.spde)
  # data stack for SPDE model
  stack.est <- inla.stack(data=list(l_counts=log(samp$counts)),
                          A=list(A.est, 1),
                          effects=list( c(mesh.index0, list(Intercept=1)),
                                        list(samp[,c("elev","trend")]) ),
                          tag='est')

  stack.pred <- inla.stack(data=list(l_counts=NA),
                           A=list(A.pred, 1),
                           effects=list( c(mesh.index0, list(Intercept=1)),
                                         list(pred[,c("elev","trend")]) ),
                           tag='pred')
  # create unified stack with tags
  stack.all <- inla.stack(stack.est, stack.pred)
  # fit model
  fit <- inla(form, 
              family="gaussian",
              data=inla.stack.data(stack.all),
              control.predictor=list(A=inla.stack.A(stack.all), compute=T, link=1),
              control.compute=c.c)
  #
    summary(fit)
  # get IDs of the test set
  idx.pred <- inla.stack.index(stack.all, "pred")$data
  
  # extract fitted values
  predvals <- fit$summary.fitted.values[idx.pred, c("0.5quant")]
  plot(pred$counts, exp(predvals),
       xlim=c(0,max(pred$counts, exp(predvals))),
       ylim=c(0,max(pred$counts, exp(predvals))))
  abline(a=0, b=1, col="red")

  plot(log(pred$counts), predvals,
     xlim=c(0,max(log(pred$counts), predvals)),
     ylim=c(0,max(log(pred$counts), predvals)))
  abline(a=0, b=1, col="red")
  
  or <- rasterFromXYZ(pred[,c("x","y","counts")], res=c(1,1))
  pr <- rasterFromXYZ(cbind(pred[,c("x","y")], exp(predvals)), res=c(1,1))

  par(mfrow=c(1,2))
  plot(or); points(samp[,c("x","y")], pch=16)
  plot(pr)
  par(mfrow=c(1,1))
  
  # extract posterior samples from fitted model
  ps <- inla.posterior.sample(10000, fit, seed=1234)
  # predict values include transformation to response + error estimate
  # psam <- sapply(ps, function(x){
  #   s.d <- sqrt(1/x$hyperpar[1])
  #   mu <- x$latent[idx.pred, 1]
  #   err <- rnorm(length(idx.pred), 0, s.d)
  #   return(exp(mu+err))
  #   ## in the GAG original scale
  # })
  # # calculate quantiles of the samples for each prediction location
  # q.sam <- t(apply(psam, 1, quantile, c(.025, 0.05, 0.5, 0.95, .975)))
  
  ps.l <- lapply(ps, function(x){ exp(x$latent[idx.pred, 1]) })
  ps.m <- do.call(rbind, ps.l)
  q.sam <- t(apply(ps.m, 2, quantile, c(0.25, 0.05, 0.5, 0.95, 0.975)))
  
  pr.s <- rasterFromXYZ(cbind(pred[,c("x","y")], q.sam[,3]), res=c(1,1))
  lci <- rasterFromXYZ(cbind(pred[,c("x","y")], q.sam[,1]), res=c(1,1))
  uci <- rasterFromXYZ(cbind(pred[,c("x","y")], q.sam[,5]), res=c(1,1))
  interval <- rasterFromXYZ(cbind(pred[,c("x","y")], (q.sam[,5]-q.sam[,1])), res=c(1,1))
  
  par(mfrow=c(1,5))
  plot(or); points(samp[,c("x","y")], pch=16)
  plot(pr.s, main="Median")
  plot(lci, main="Lower CI")
  plot(uci, main="Upper CI")
  plot(interval, main="Range")
  par(mfrow=c(1,1))
  
  plot(pred$counts, q.sam[,3],
       xlim=c(0,max(pred$counts, q.sam[,3])),
       ylim=c(0,max(pred$counts, q.sam[,3])))
  abline(a=0, b=1, col="red")
  
  # total
  p.tot <- apply(psam, 2, sum)
  hist(p.tot); summary(p.tot)
}

# fit$formula <- form
# pi <- inlabru:::predict.inla(fit, pred, ~ -1 + Intercept + elev + trend)
