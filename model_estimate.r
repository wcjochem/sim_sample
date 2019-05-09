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
  form <- pop ~ -1 + Intercept + cov + f(sett, model="iid", values=1:5) + f(field, model=spde)
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

  # stack.pred_in <- inla.stack(data=list(pop=NA),
  #                             A=list(A.pred_in, 1),
  #                             effects=list( c(mesh.index0, list(Intercept=1)),
  #                                           list(pred_in[,c("cov","sett")]) ),
  #                             tag='pred_in')
  # 
  # stack.pred_out <- inla.stack(data=list(pop=NA),
  #                              A=list(A.pred_out, 1),
  #                              effects=list( c(mesh.index0, list(Intercept=1)),
  #                                            list(pred_out[,c("cov","sett")]) ),
  #                              tag='pred_out')
    
  # create unified stack with tags
  # stack.all <- inla.stack(stack.est, stack.pred)
  # fit model
  fit <- tryCatch({
      inla(form, 
           family="poisson",
           data=inla.stack.data(stack.est),
           control.predictor=list(A=inla.stack.A(stack.est), compute=T, link=1),
           control.compute=c.c, verbose=T)
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
    nsamp <- 1e4
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
    # within sample area
    linpred <- as.matrix(A.pred_in %*% xSpace + xSett[pred_in$sett,] + as.matrix(cbind(1, pred_in[,c("cov")])) %*% xX)
    pred_N_in <- t(apply(linpred, 1, FUN=function(x){ quantile(rpois(n=nsamp, lambda=exp(x)), probs=c(0.025,0.5,0.975)) }))
    # outside sample area
    linpred <- as.matrix(A.pred_out %*% xSpace + xSett[pred_out$sett,] + as.matrix(cbind(1, pred_out[,c("cov")])) %*% xX)
    pred_N_out <- t(apply(linpred, 1, FUN=function(x){ quantile(rpois(n=nsamp, lambda=exp(x)), probs=c(0.025,0.5,0.975)) }))
  }
  ## return
  return(list("predN_in"=pred_N_in, "predN_out"=pred_N_out))
  
}