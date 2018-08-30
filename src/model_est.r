#
# Code for simplified model-based estimates
# Used with sim_sample_test.r 
#
# Chris Jochem
# 05 June 2018
#

# require(inlabru) 
require(INLA)
require(gbm)
require(party)
require(MASS)
#require(dismo)

## replicate ORNL confidence interval estimation ##
# Crude Monte Carlo from SRS -- shown in DRC work #
mc_ci <- function(samp=NULL, total_area=NULL){
  sampsize <- nrow(samp) # samples; total_area is the length of the domain

  total_est <- mean(samp$counts * total_area) # average predicted total
  sq_err <- ((samp$counts * total_area) - total_est)^2 # squared errors
  sum_sqerr <- sum(sq_err) # sum of squared errors
  # estimate the variance
  var_est <- sum_sqerr / (sampsize * (sampsize - 1))
  t_95 <- qt(0.975, sampsize - 1) # approximate using t distribution
  margins <- t_95 * var_est^.5 # upper/lower CI
  return(c(total_est - margins, total_est + margins))
}

# lognormal method from Weber et al. paper #
lognorm_ci <- function(samp=NULL, total_area=NULL){
  sampsize <- nrow(samp) # samples; total_area is the length of the domain
  
  total_est <- vector("numeric",length=10000) # store 10k realisations
  fit_param <- fitdistr(samp$counts, "lognormal") # lognormal dist fit to the samples
  for(i in 1:10000){ # draw 10k different approximate samples with random alternation around lognorm
    m_est <- rnorm(1,mean = fit_param$estimate["meanlog"], fit_param$sd["meanlog"]) # mean
    sd_est <- rnorm(1,mean = fit_param$estimate["sdlog"], fit_param$sd["sdlog"]) # sd
    total_est[i] <- sum(rlnorm(total_area, m_est, sd_est)) # draw density value per pixel from the lognormal
  }
  return(quantile(total_est, probs=c(0.025, 0.975)))
}

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
# model-based geostatistical function
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
  # stack.est <- inla.stack(data=list(l_counts=log(samp$counts)),
  #                         A=list(A.est, 1),
  #                         effects=list( c(mesh.index0, list(Intercept=1)),
  #                                       list(samp[,c("elev","trend")]) ),
  #                         tag='est')
  
  stack.est <- inla.stack(data=list(l_counts=samp$counts),
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
  # fit <- inla(form, 
  #             family="gaussian",
  #             data=inla.stack.data(stack.all),
  #             control.predictor=list(A=inla.stack.A(stack.all), compute=T, link=1),
  #             control.compute=c.c)
  #
  fit <- tryCatch({
      inla(form, 
           family="nbinomial",
           data=inla.stack.data(stack.all),
           control.predictor=list(A=inla.stack.A(stack.all), compute=T, link=1),
           control.compute=c.c)
    },
    error=function(e){
      NULL
    }
  )
  # check return
  if(is.null(fit)){
    predvals <- NA
  } else{
    # summary(fit)
    # get IDs of the test set
    idx.pred <- inla.stack.index(stack.all, "pred")$data
    # extract fitted values - posterior median of response
    # predvals <- exp(fit$summary.fitted.values[idx.pred, c("0.5quant")])
    predvals <- fit$summary.fitted.values[idx.pred, c("0.5quant")]
  }
  ## return
  return(list("predvals"=predvals, "fittedmod"=fit))
  
  #####
  # plot(pred$counts, exp(predvals),
  #      xlim=c(0,max(pred$counts, exp(predvals))),
  #      ylim=c(0,max(pred$counts, exp(predvals))))
  # abline(a=0, b=1, col="red")
  # 
  # plot(log(pred$counts), predvals,
  #    xlim=c(0,max(log(pred$counts), predvals)),
  #    ylim=c(0,max(log(pred$counts), predvals)))
  # abline(a=0, b=1, col="red")
  # 
  # plot(pred$counts, predvals,
  #      xlim=c(0,max(pred$counts, predvals)),
  #      ylim=c(0,max(pred$counts, predvals)))
  # abline(a=0, b=1, col="red")
  # cor(pred$counts, predvals)^2
  # 
  # or <- rasterFromXYZ(pred[,c("x","y","counts")], res=c(1,1))
  # pr <- rasterFromXYZ(cbind(pred[,c("x","y")], exp(predvals)), res=c(1,1))
  # pr <- rasterFromXYZ(cbind(pred[,c("x","y")], predvals), res=c(1,1))
  # 
  # par(mfrow=c(1,2))
  # plot(or); points(samp[,c("x","y")], pch=16)
  # plot(pr)
  # par(mfrow=c(1,1))
  # 
  # # extract posterior samples from fitted model
  # ps <- inla.posterior.sample(10000, fit, seed=1234)
  # # back-transform
  # ps.l <- lapply(ps, function(x){ exp(x$latent[idx.pred, 1]) })
  # ps.m <- do.call(rbind, ps.l)
  # q.sam <- t(apply(ps.m, 2, quantile, c(0.25, 0.05, 0.5, 0.95, 0.975)))
  # 
  # pr.s <- rasterFromXYZ(cbind(pred[,c("x","y")], q.sam[,3]), res=c(1,1))
  # lci <- rasterFromXYZ(cbind(pred[,c("x","y")], q.sam[,1]), res=c(1,1))
  # uci <- rasterFromXYZ(cbind(pred[,c("x","y")], q.sam[,5]), res=c(1,1))
  # interval <- rasterFromXYZ(cbind(pred[,c("x","y")], (q.sam[,5]-q.sam[,1])), res=c(1,1))
  # 
  # par(mfrow=c(1,5))
  # plot(or); points(samp[,c("x","y")], pch=16)
  # plot(pr.s, main="Median")
  # plot(lci, main="Lower CI")
  # plot(uci, main="Upper CI")
  # plot(interval, main="Range")
  # par(mfrow=c(1,1))
  # 
  # plot(pred$counts, q.sam[,3],
  #      xlim=c(0,max(pred$counts, q.sam[,3])),
  #      ylim=c(0,max(pred$counts, q.sam[,3])))
  # abline(a=0, b=1, col="red")
  # 
  # # total
  # p.tot <- apply(psam, 2, sum)
  # hist(p.tot); summary(p.tot)
}

# boosted regression tree
brt <- function(samp, pred=NULL){
  # log transform outcome
  samp$l_counts <- log(samp$counts)
  # 
  # fit <- gbm.step(data=samp,
  #                 gbm.x=c("elev","trend","x","y"),
  #                 gbm.y="l_counts",
  #                 family="gaussian",
  #                 tree.complexity=4,
  #                 learning.rate=0.001,
  #                 bag.fraction=0.5)
  
  # fit model on sample
  fit <- tryCatch({
      gbm(l_counts ~ elev + trend,
          distribution="gaussian",
          data=samp,
          n.trees=2000,
          interaction.depth=4,
          shrinkage=0.001,
          bag.fraction=0.5)
  },
  error=function(x){
    NULL
  })
  # check return
  if(is.null(fit)){
    predvals <- NA
  } else{
    # make predictions into the full domain
    predvals <- predict.gbm(fit,
                            newdata=pred,
                            n.trees=2000,
                            type="link")
    predvals <- exp(predvals) # back-transform to counts
  }
  ## return
  return(list("predvals"=predvals, "fittedmod"=fit))
}

# random forest model
rf <- function(samp, pred=NULL){
  # transform outcome
  samp$l_counts <- log(samp$counts)
  
  # fit model on sample
  fit <- tryCatch({
    cforest(l_counts ~ elev + trend,
            data=samp,
            controls=cforest_unbiased(mtry=2, ntree=1000))
  },
  error=function(e){
    NULL
  })
  # check return
  if(is.null(fit)){
    predvals <- NA
  } else{
    # make predictions
    predvals <- predict(fit,
                        pred,
                        OOB=TRUE,
                        type="response")
    # back transform
    predvals <- exp(as.vector(predvals))
  }
  ## return
  return(list("predvals"=predvals, "fittedmod"=fit))
}


