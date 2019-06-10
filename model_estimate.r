# 
# Estimate a variety of simplified population models
# model_estimate.r
#
# Helper functions for model-fitting
#
# Author: Chris Jochem (W.C.Jochem@soton.ac.uk)
#


# Model-based geostatistics
mbg <- function(samp=NULL, pred_in=NULL, pred_out=NULL, mesh=NULL){
  # spde - spatial prior
  spde <- inla.spde2.pcmatern(mesh, prior.range=c(10, .9), prior.sigma=c(.5, .5))
  # components
  form <- pop ~ -1 + Intercept + cov + f(sett, model="iid") + f(field, model=spde)
  c.c <- list(cpo=TRUE, dic=TRUE, waic=TRUE, config=TRUE)
  # set up model
  A.est <- inla.spde.make.A(mesh=mesh,
                            loc=data.matrix(samp[,c("x","y")]))
  
  A.pred_in <- inla.spde.make.A(mesh=mesh,
                                loc=data.matrix(pred_in[,c("x","y")]))
  
  A.pred_out <- inla.spde.make.A(mesh=mesh,
                                 loc=data.matrix(pred_out[,c("x","y")]))
  # index to the mesh
  mesh.index0 <- inla.spde.make.index(name="field", n.spde=spde$n.spde)
  # data stack for model
  stack.est <- inla.stack(data=list(pop=samp$pop),
                          A=list(A.est, 1),
                          effects=list( c(mesh.index0, list(Intercept=1)),
                                        list(samp[,c("cov","sett")]) ),
                          tag='est')
  # fit model
  fit <- tryCatch({
      inla(form, 
           family="poisson",
           data=inla.stack.data(stack.est),
           control.predictor=list(A=inla.stack.A(stack.est), compute=T, link=1),
           control.compute=c.c 
           # verbose=T
          )
    },
    error=function(e){
      NULL
    }
  )
  # check return
  if(is.null(fit)){
    predvals <- NA
  } else{
    # predictions
    nsamp <- 1e3
    # draw samples from the posterior
    ps <- inla.posterior.sample(n=nsamp, fit)
    # get indices to the effects
    contents <- fit$misc$configs$contents
    
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
    linpred <- as.matrix(A.est %*% xSpace + xSett[samp$sett,] + as.matrix(cbind(1, samp[,c("cov")])) %*% xX)
    pred_N_s <- t(apply(linpred, 1, FUN=function(x){ quantile(rpois(n=nsamp, lambda=exp(x)), probs=c(0.025,0.5,0.975)) }))
    # within sample area
    linpred <- as.matrix(A.pred_in %*% xSpace + xSett[pred_in$sett,] + as.matrix(cbind(1, pred_in[,c("cov")])) %*% xX)
    pred_N_in <- t(apply(linpred, 1, FUN=function(x){ quantile(rpois(n=nsamp, lambda=exp(x)), probs=c(0.025,0.5,0.975)) }))
    # outside sample area
    linpred <- as.matrix(A.pred_out %*% xSpace + xSett[pred_out$sett,] + as.matrix(cbind(1, pred_out[,c("cov")])) %*% xX)
    pred_N_out <- t(apply(linpred, 1, FUN=function(x){ quantile(rpois(n=nsamp, lambda=exp(x)), probs=c(0.025,0.5,0.975)) }))
  }
  ## return
  return(list("predN_s"=pred_N_s,"predN_in"=pred_N_in, "predN_out"=pred_N_out))
  
}


# Weighted Model-based geostatistics
wmbg <- function(samp=NULL, pred_in=NULL, pred_out=NULL, mesh=NULL, wgts=NULL){
  # spde - spatial prior
  spde <- inla.spde2.pcmatern(mesh, prior.range=c(10, .9), prior.sigma=c(.5, .5))
  # components
  form <- pop ~ -1 + Intercept + cov + f(sett, model="iid") + f(field, model=spde)
  c.c <- list(cpo=TRUE, dic=TRUE, waic=TRUE, config=TRUE)
  # set up model
  A.est <- inla.spde.make.A(mesh=mesh,
                            loc=data.matrix(samp[,c("x","y")]))
  
  A.pred_in <- inla.spde.make.A(mesh=mesh,
                                loc=data.matrix(pred_in[,c("x","y")]))
  
  A.pred_out <- inla.spde.make.A(mesh=mesh,
                                 loc=data.matrix(pred_out[,c("x","y")]))
  # index to the mesh
  mesh.index0 <- inla.spde.make.index(name="field", n.spde=spde$n.spde)
  # data stack for model
  stack.est <- inla.stack(data=list(pop=samp$pop,
                                    invwts=((1/wgts)/sum(1/wgts))*length(wgts)), #1/wgts ;((1/wgts) / sum(1/wgts))
                          A=list(A.est, 1),
                          effects=list( c(mesh.index0, list(Intercept=1)),
                                        list(samp[,c("cov","sett")]) ),
                          tag='est')
  # fit model
  inla.setOption("enable.inla.argument.weights", TRUE)
  fit <- tryCatch({
      inla(form, 
           family="poisson",
           data=inla.stack.data(stack.est),
           weights=inla.stack.data(stack.est)$invwts,
           control.predictor=list(A=inla.stack.A(stack.est), compute=T, link=1),
           control.compute=c.c
           #verbose=T
          )
    },
    error=function(e){
      NULL
    }
  )
  # check return
  if(is.null(fit)){
    predvals <- NA
  } else{
    # predictions
    nsamp <- 1e3
    # draw samples from the posterior
    ps <- inla.posterior.sample(n=nsamp, fit)
    # get indices to the effects
    contents <- fit$misc$configs$contents
    
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
    linpred <- as.matrix(A.est %*% xSpace + xSett[samp$sett,] + as.matrix(cbind(1, samp[,c("cov")])) %*% xX)
    pred_N_s <- t(apply(linpred, 1, FUN=function(x){ quantile(rpois(n=nsamp, lambda=exp(x)), probs=c(0.025,0.5,0.975)) }))
    # within sample area
    linpred <- as.matrix(A.pred_in %*% xSpace + xSett[pred_in$sett,] + as.matrix(cbind(1, pred_in[,c("cov")])) %*% xX)
    pred_N_in <- t(apply(linpred, 1, FUN=function(x){ quantile(rpois(n=nsamp, lambda=exp(x)), probs=c(0.025,0.5,0.975)) }))
    # outside sample area
    linpred <- as.matrix(A.pred_out %*% xSpace + xSett[pred_out$sett,] + as.matrix(cbind(1, pred_out[,c("cov")])) %*% xX)
    pred_N_out <- t(apply(linpred, 1, FUN=function(x){ quantile(rpois(n=nsamp, lambda=exp(x)), probs=c(0.025,0.5,0.975)) }))
  }
  ## return
  return(list("predN_s"=pred_N_s,"predN_in"=pred_N_in, "predN_out"=pred_N_out))
}

# preferential sample mbg
prefsamp <- function(){
  # spde - spatial prior
  spde <- inla.spde2.pcmatern(mesh, prior.range=c(10, .9), prior.sigma=c(.5, .5))
  # components
  form <- pop ~ -1 + Intercept + cov + f(sett, model="iid") + f(field, model=spde)
  c.c <- list(cpo=TRUE, dic=TRUE, waic=TRUE, config=TRUE)
  # set up model
  A.est <- inla.spde.make.A(mesh=mesh,
                            loc=data.matrix(samp[,c("x","y")]))
  
  A.pred_in <- inla.spde.make.A(mesh=mesh,
                                loc=data.matrix(pred_in[,c("x","y")]))
  
  A.pred_out <- inla.spde.make.A(mesh=mesh,
                                 loc=data.matrix(pred_out[,c("x","y")]))
  # index to the mesh
  mesh.index0 <- inla.spde.make.index(name="field", n.spde=spde$n.spde)
  # data stack for model
  # stack.est <- inla.stack(data=list(pop=samp$pop,
  #                                   invwts=1/wgts), #((1/wgts) / sum(1/wgts))
  #                         A=list(A.est, 1),
  #                         effects=list( c(mesh.index0, list(Intercept=1)),
  #                                       list(samp[,c("cov","sett")]) ),
  #                         tag='est')
  
  stk.y <- inla.stack(data=list(y=cbind(pop=samp$pop, NA), e=rep(0,nrow(samp))),
                      A=list(A.est, 1),
                      effects=list(i=1:mesh$n, Intercept=rep(1,nrow(samp))),
                      tag='resp')
  
  stk.pp <- inla.stack(data=list(y=cbind(NA, y.pp), e=e.pp))
  

stk2.pp <- inla.stack(data = list(y = cbind(NA, y.pp), e = e.pp), 
  A = list(A.pp, 1),
  effects = list(j = 1:nv, b0.pp = rep(1, nv + n)),
  tag = 'pp2')

j.stk <- inla.stack(stk2.y, stk2.pp)
  
  
}
