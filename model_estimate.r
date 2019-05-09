# 
# Estimate a variety of simplified population models
# model_estimate.r
#
# Helper functions for model-fitting
#
# Author: Chris Jochem (W.C.Jochem@soton.ac.uk)
#

# Model-based geostatistics
mbg <- function(samp=NULL, pred=NULL, mesh=NULL){
  # spde - spatial prior
  spde <- inla.spde2.pcmatern(mesh, prior.range=c(10, .9), prior.sigma=c(.5, .5))
  # components
  form <- pop ~ -1 + Intercept + cov + f(field, model=spde)
  c.c <- list(cpo=TRUE, dic=TRUE, waic=TRUE, config=TRUE)
  # set up model
  A.est <- inla.spde.make.A(mesh=mesh,
                            loc=data.matrix(samp[,c("x","y")]))
  
  A.pred <- inla.spde.make.A(mesh=mesh,
                             loc=data.matrix(pred[,c("x","y")]))
  # index to the mesh
  mesh.index0 <- inla.spde.make.index(name="field", n.spde=spde$n.spde)
  # data stack for model
  stack.est <- inla.stack(data=list(pop=samp$pop),
                          A=list(A.est, 1),
                          effects=list( c(mesh.index0, list(Intercept=1)),
                                        list(samp[,c("cov","sett")]) ),
                          tag='est')

  stack.pred <- inla.stack(data=list(pop=NA),
                           A=list(A.pred, 1),
                           effects=list( c(mesh.index0, list(Intercept=1)),
                                         list(pred[,c("cov","sett")]) ),
                           tag='pred')
  # create unified stack with tags
  stack.all <- inla.stack(stack.est, stack.pred)
  # fit model
  fit <- tryCatch({
      inla(form, 
           family="poisson",
           data=inla.stack.data(stack.all),
           control.predictor=list(A=inla.stack.A(stack.all), compute=T, link=1),
           control.compute=c.c, verbose=F)
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
  
}