##  A fac-ipm simulator based exactly on the 'Spirit' version of the model: demographic componenets enter as informed priors.
##
##
##  05 jul 2022
####################################################

## Loads a function, that's all.

gen.noCovs <- function(
  n.sites = c(100,100),                 ## br, nb
  n.yr = 10,
  n.sim = 100,
  phi.mn = c(0.9, 0.75, 0.75, 0.8),    ## mean survival: br, aut, nb, spr
  phi.sd = rep(0,4),  #c(0.2, 0.3, 0.2, 0.3),      ## annual variation in survival: br, aut, nb, spr
  rho.mn = 1.1,   ## 0.3860718,                  ## natural scale: expected value from the `br.dmgr.sum` in `dat` --> N.b. this yields 'per pair' reproduction.  Per capita rho is this divided by two!
  rho.sd = 0, #0.5,                        ## reproduction: log scale
  eps.mn = c(0.2, 0.9),                ## emigration: br, nb
  iota.mn = c(0.65, 0.52),                 ## immigration: br, nb
  psi.mn = NULL, #c(0.7, 0.7),                ## occupancy (constant across years): br, nb; or NULL if you don't want any zero inflation.
  lambda.mn = c(5, 1),              ## initial abundance: br, nb
  p = c(0.7, 0.7)                      ## detection probability: br, nb
) {
 
  
  require(mvtnorm)

  br.dim <- c(n.sites[1], n.yr, n.sim)
  nb.dim <- c(n.sites[2], n.yr, n.sim)


### covariate values are the same for all simulations:
#if (!is.null(br.sigma)) {
#  br.covs <- array(rmvnorm(n = n.sites[1] * n.yr, sigma = br.sigma),
#                 dim = c(n.sites[1], n.yr, nrow(br.sigma)))
#}

#if (!is.null(nb.sigma)) {
#  nb.covs <- array(rmvnorm(n = n.sites[2] * n.yr, sigma = nb.sigma),
#                 dim = c(n.sites[2], n.yr, nrow(nb.sigma)))
#}


  phi.int <- qlogis(phi.mn)


#stopifnot(prod(phi.mn[1:3]) > phi.ann.mn)

## at first, no covariance:
  phi <- aperm(array(plogis(rnorm(n = 4*n.yr*n.sim, mean = phi.int, sd = phi.sd)),
              dim = c(4, n.yr, n.sim),
              ), perm = c(2,1,3))

  emp.phi.mn <- apply(phi, MAR = 2:3, FUN = mean)
  emp.phi.sd <- apply(qlogis(phi), MAR = 2:3, FUN = sd)

  phi.ann <- apply(phi, MAR = c(1,3), FUN = prod)


  rho.int <- log(rho.mn)

  rho <- matrix(exp(rnorm(n = n.yr*n.sim, mean = log(rho.mn), sd = rho.sd)),
              nrow = n.yr, ncol = n.sim)

  emp.rho.mn <- apply(rho, MAR = 2, FUN = mean)
  
## constant expected value:
  eps.int <- qlogis(eps.mn)
  eps <- plogis(eps.int)

## constant expected value:
  iota.int <- log(iota.mn)
  iota <- exp(iota.int)


##  Occupancy
  if (!is.null(psi.mn)) {
    psi.int <- qlogis(psi.mn)
    psi <- plogis(psi.int)
    z.br <- matrix(rbinom(n = n.sites[1]*n.sim, prob = psi[1], size = 1), nrow = n.sites[1], ncol = n.sim, dimnames = list(site = NULL, sim = NULL))
    z.nb <- matrix(rbinom(n = n.sites[2]*n.sim, prob = psi[2], size = 1), nrow = n.sites[1], ncol = n.sim, dimnames = list(site = NULL, sim = NULL))
  } else {
    psi.int <- psi <- NA
    z.br <- matrix(1, nrow = n.sites[1], ncol = n.sim, dimnames = list(site = NULL, sim = NULL))
    z.nb <- matrix(1, nrow = n.sites[2], ncol = n.sim, dimnames = list(site = NULL, sim = NULL))
  }
  


##  Initial abundance
  emp.eps.br <- b.br <- S.br <- G.br <- I.br <- N.br <- array(NA, dim = br.dim, dimnames = list(site = NULL, yr = NULL, sim = NULL))
  emp.eps.nb <- b.nb <- S.nb <- G.nb <- I.nb <- N.nb <- array(NA, dim = nb.dim, dimnames = list(site = NULL, yr = NULL, sim = NULL))

  lambda.int <- log(lambda.mn)
  lambda <- exp(lambda.int)

  N.br[,1,] <- rpois(n = n.sites[1]*n.sim, lambda = z.br * lambda[1])
  N.nb[,1,] <- rpois(n = n.sites[2]*n.sim, lambda = z.nb * lambda[2])

  gamma.nb <- phi.nb <- gamma.br <- phi.br <- matrix(NA, nrow = n.yr, ncol = n.sim, dimnames = list(yr = NULL, sim = NULL))

  for (yr in 2:n.yr) {

##  Dynamics: breeding_____________________________________________________________________________________________
    phi.br[yr,] <- apply(phi[yr-1,,], MAR = 2, FUN = prod)
    gamma.br[yr,] <- rho[yr-1,] * apply(phi[yr-1,2:4,], MAR = 2, FUN = prod) 
## effective number of breeders!
    b.br[,yr,] <- N.br[,yr-1,] * (1 - eps[1]) + iota[1] / phi.br[yr,]
    
    S.br[,yr,] <- rbinom(n = n.sites[1]*n.sim, size = N.br[,yr-1,], prob = rep(phi.br[yr,]*(1 - eps[1]), each = n.sites[1]))
    I.br[,yr,] <- rpois(n = n.sites[1]*n.sim, lambda = z.br * iota[1])

    G.br[,yr,] <- rpois(n = n.sites[1]*n.sim, lambda = rep(gamma.br[yr,], each = n.sites[1]) * b.br[,yr,])

    N.br[,yr,] <- S.br[,yr,] + G.br[,yr,] + I.br[,yr,]

    emp.eps.br[,yr,] <- 1 - (S.br[,yr,]/(max(N.br[,yr-1,], 0.001)*phi.br[yr,]))
## non-breeding:_____________________________________________________________________________________________
    phi.nb[yr,] <- apply(phi[yr-1,3:4,], MAR = 2, FUN = prod) * apply(phi[yr,1:2,], MAR = 2, FUN = prod)
    gamma.nb[yr,] <- rho[yr,]/2 * phi[yr,2,]

## effective number of breeders!
    b.nb[,yr,] <- N.nb[,yr-1,] * apply(phi[yr-1,3:4,], MAR = 2, FUN = prod) * (1 - eps[2]) + (iota[2] / apply(phi[yr,1:2,], MAR = 2, FUN = prod))

    S.nb[,yr,] <- rbinom(n = n.sites[2], size = N.nb[,yr-1,], prob = rep(phi.nb[yr]*(1 - eps[2]), each = n.sites[2]))
    I.nb[,yr,] <- rpois(n = n.sites[2], lambda = z.nb * iota[2])
    G.nb[,yr,] <- rpois(n = n.sites[2], lambda = rep(gamma.nb[yr,], each = n.sites[2]) * b.nb[,yr,])

    N.nb[,yr,] <- S.nb[,yr,] + G.nb[,yr,] + I.nb[,yr,]

    emp.eps.nb[,yr,] <- 1 - (S.nb[,yr,]/(max(N.nb[,yr-1,], 0.001)*phi.nb[yr,]))

  }

  obs.br <- array(rbinom(n = prod(dim(N.br)), size = N.br, prob = p[1]), dim = dim(N.br))
  obs.nb <- array(rbinom(n = prod(dim(N.nb)), size = N.nb, prob = p[2]), dim = dim(N.nb))

  grwth.br <- {tmp <- apply(N.br, MAR = 2:3, FUN = mean); apply(tmp[2:nrow(tmp),]/tmp[1:(nrow(tmp)-1),], MAR = 2, FUN = function(v) { prod(v)^(1/length(v)) } )}
  grwth.nb <- {tmp <- apply(N.nb, MAR = 2:3, FUN = mean); apply(tmp[2:nrow(tmp),]/tmp[1:(nrow(tmp)-1),], MAR = 2, FUN = function(v) { prod(v)^(1/length(v)) } )}  
  
  N.br.mn <- apply(N.br, MAR = 2:3, FUN = mean)
  S.br.mn <- apply(S.br, MAR = 2:3, FUN = mean, na.rm=T)


  N.nb.mn <- apply(N.nb, MAR = 2:3, FUN = mean)
  S.nb.mn <- apply(S.nb, MAR = 2:3, FUN = mean, na.rm=T)
  
  emp.eps <- c(mean(1 - (S.br.mn[2:nrow(S.br.mn),]/(N.br.mn[1:(nrow(N.br.mn)-1),]*phi.br[2:nrow(phi.br),]))),
                   mean(1 - (S.nb.mn[2:nrow(S.nb.mn),]/(N.nb.mn[1:(nrow(N.nb.mn)-1),]*phi.nb[2:nrow(phi.nb),]))))

  
  return(
  list(dat = mget(c("obs.br", "obs.nb", "n.yr", "n.sites", "p")),
    true = mget(c("psi.int", "lambda.int", "phi.int", "phi.sd", "eps.int", "iota.int", "rho.int", "rho.sd", "psi", "phi", "rho")),
    emp = mget(c(phi.mn = "emp.phi.mn", phi.sd = "emp.phi.sd", rho.mn = "emp.rho.mn", eps.mn = "emp.eps")),
    grwth = mget(c("grwth.br", "grwth.nb")))
  )

}  ## END gen.noCovs definition



