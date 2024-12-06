##  Standard functions to run JAGS models in parallel.
##
##
##  07 feb 2022
#############################################

######################################################################
##  BEGIN DEFINE `jags.fun`, which runs a single chain of a JAGS model.
######################################################################
jags.fun <- function (l) {

  require(parallel)
  require(rjags)

  if (!all(c("glm", "lecuyer") %in% rjags::list.modules())) {
    load.module("glm")
    load.module("lecuyer")
  }
    
  init.fun <- function() {
    c(l$inits,
      list(
      .RNG.state = l$.RNG.state,
      .RNG.name = l$.RNG.name))
  }
  
  jm <- jags.model(file=l$file.nm, inits=init.fun(), data=l$dat, n.adapt = l$adapt, n.chains=1) 
  update(jm, n.iter=l$iter.burn)
  js <- coda.samples(jm, var=l$var, n.iter=l$iter.samp, thin = l$thin)
  want <- c("js", "jm")
  return(mget(want))
}
######################################################################
##  END DEFINE `jags.fun`
######################################################################




#################################################################################
##  BEGIN DEFINE `jags.cont`, which picks up running a previously run JAGS model.
#################################################################################
jags.cont <- function (jm, n.update, n.monitor, vars) {
  jm$recompile()
  rjags:::update.jags(jm, n.iter=n.update)
  js <- coda.samples(jm, var=vars, n.iter=n.monitor)
  want <- c("js", "jm")
  return(mget(want))
}  
#################################################################################
##  END DEFINE `jags.cont`
#################################################################################

