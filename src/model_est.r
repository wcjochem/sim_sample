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
# make mesh
makemesh <- function(locn=NULL, bound=NULL){
  mesh <- inla.mesh.2d(loc=locn,
                       boundary=bound,
                       max.edge=c(2, 5),
                       cutoff=1,
                       offset=c(5, 10))

  return(mesh)
}
