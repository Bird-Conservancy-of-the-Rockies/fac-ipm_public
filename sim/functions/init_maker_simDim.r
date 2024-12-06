##  Function to make initial values for fitting JAGS models to fac-ipm simulations.
##  Here, resulting objects will have a third dimension: n.sim.
##
##  Jul 2022, BLN
###################################

### 1 --
 
arr.expander <- function (a, new.dim) {    ## copies the old array 'new.dim' (or 'prod(new.dim)', if a vector)number of times, to fill in a new array of one greater dimension, with that dimension = new.dim. 
  
  if (is.vector(a)) {
    old.dim <- length(a)
  } else {
    old.dim <- dim(a)
  }
  
  array(rep(c(a), new.dim), dim = c(old.dim, new.dim))
}

#########################################################################

### 2 -- 

init.maker <- function (
  dat.lst,           ##  the data list containing the `br.obs` and `nb.obs` arrays: these will have a third dimension, n.sim.
  simDim = FALSE     ##  should the created objects be copied out into larger arrays that have their final dimension equal to 'n.sim' inside of 'dat.lst'
) {

  with (dat.lst, {
    n.sim <- dim(obs.br)[3]
    est.N.br.max <- round(max(obs.br)/p[1])
    est.N.nb.max <- round(max(obs.nb)/p[2])


    out.lst <<- list(
      N.br = cbind(rep(est.N.br.max, n.sites[1]), matrix(NA, nrow = n.sites[1], ncol = n.yr-1)),
      S.br = cbind(NA, matrix(est.N.br.max, nrow = n.sites[1], ncol = n.yr-1)),
      G.br = cbind(NA, matrix(est.N.br.max, nrow = n.sites[1], ncol = n.yr-1)),
#      E.br = cbind(NA, matrix(0, nrow = n.sites[1], ncol = n.yr-1)),
      I.br = cbind(NA, matrix(0, nrow = n.sites[1], ncol = n.yr-1)),
      z.br = rep(1, n.sites[1]),

      N.nb = cbind(rep(est.N.nb.max, n.sites[2]), matrix(NA, nrow = n.sites[2], ncol = n.yr-1)),
      S.nb = cbind(NA, matrix(est.N.nb.max, nrow = n.sites[2], ncol = n.yr-1)),
      G.nb = cbind(NA, matrix(est.N.nb.max, nrow = n.sites[2], ncol = n.yr-1)),
#      E.nb = cbind(NA, matrix(0, nrow = n.sites[2], ncol = n.yr-1)),
      I.nb = cbind(NA, matrix(0, nrow = n.sites[2], ncol = n.yr-1)),
      z.nb = rep(1, n.sites[2])
    )
    
  if (simDim) {
    out.lst <<- rapply(out.lst, f = arr.expander, new.dim = n.sim, how = "replace") 
  }
  
  })  ## END with(dat.lst, ...)
 
  return(out.lst) 
}
