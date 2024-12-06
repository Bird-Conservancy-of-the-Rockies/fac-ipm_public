##  Parallel Dail-Madsen models, with Kuo-Mallick variable selection.
##  A one-species version!
##
##  25 may 2022
###################################

##  Q: How do Year labels compare with Transition labels?
##  A: They match!  That is, transition from Year 1 to Year 2 is called phi.xx.ann[1] ! 
##  A2:  And because the years have already been lined up properly between breeding and non-breeding counts, the survival from breeding Year 1 to Year 2 is called phi.br.ann[1], and the survival from non-breeding Year 1 to Year 2 is called phi.nb.ann[1] ! 

data {

## Standard log-linear and logit-linear priors for normal distribution precision:  
  log.tau <- 0.05
  logit.tau <- 0.37   ## Jeffreys prior

  small <- 0.001

  K.br <- 400

#  no.hyper <- c(2,4)
#  yes.hyper <- c(1,3)
#  
#  phi.hyper.mn <- c(0.9, -99, 0.85, -99)
#  phi.hyper.sd <- c(0.2, -99, 0.2, -99)

}  ## END data

model {

for (si in 1:n.sim) {

####################################################################
## Priors:
br.phi.mn[si] ~ dbeta(9,1)  ## beta is further from symmetrical, but realistic way to apply a prior. 
phi.int[1,si] <- logit(br.phi.mn[si])  #~ dnorm(phi.int.true[1], (phi.int.true[1]/3)^(-2))   ## CV = 0.33

#nb.phi.mn[si] ~ dbeta(7,3)  ## beta is further from symmetrical, but realistic way to apply a prior. 
#phi.int[3,si] <- logit(nb.phi.mn[si])  #~ dnorm(phi.int.true[1], (phi.int.true[1]/3)^(-2))   ## CV = 0.33

#phi.int[2,si] ~ dnorm(0,logit.tau)
#phi.int[4,si] ~ dnorm(0,logit.tau)

  for (q in 2:4) {  ## Loop over seasons: order is 1) summer, 2) fall, 3) winter, 4) spring. 
    phi.int[q,si] ~ dnorm(0,logit.tau)
#    phi.sd[q,si] ~ dunif(0,1.64399)  ## 1.64399 = logit.tau^(-0.5)  
  }

rho.int[si] ~ dnorm(0.095, 0.08^(-2))
#  rho.sd[si] ~ dunif(0,log.tau^(-0.5))  ## 1.64399 = logit.tau^(-0.5)  

for (yr in 1:n.yr) {  ## Loop over years (add one at the end of the series of counts)

for (q in 1:4) {  ## Loop over seasons: order is 1) summer, 2) fall, 3) winter, 4) spring. 
#    phi.yr[yr,q,si] ~ dnorm(0, phi.sd[q,si]^(-2))
    logit(phi[yr,q,si]) <- phi.int[q,si] #+ phi.yr[yr,q,si]
}

  phi.ann[yr,si] <- prod(phi[yr,,si])

## pre-breeding census: beginning of summer.
  grwth.pop[yr,si] <- phi.ann[yr,si] + (prod(phi[yr,2:4,si]) * rho[yr,si]/2)  ## the '2' translates per male to per capita reproduction!

#  rho.yr[yr,si] ~ dnorm(0, rho.sd[si]^(-2))
  log(rho[yr,si]) <- rho.int[si]# + rho.yr[yr,si]
  
}  ## END loop over years


############################################################################
##  Linear predictor intercepts: it's just easier to specify these by name here.

for (q in 1:2) {  ## here, 1 = br, 2 = nb
#  psi.int[q,si] ~ dnorm(0,logit.tau)
  eps.int[q,si] ~ dnorm(0,logit.tau) #T(,0)   ## ADDED 20 Jun 2022
  lamb.int[q,si] ~ dnorm(0,log.tau)
  iota.int[q,si] ~ dnorm(0,log.tau)
}


for (ce in 1:n.sites[1]) {
  
#  log(lambda.br[ce]) <- lamb.br.int #+
#    beta.lamb.br[1] * br.x.arr[ce,1,6] +  #1 prcp
#    alpha.lamb.br[br.ce.strat[ce]]  ## random stratum effect
#  
#  z.br[ce,si] ~ dbern(ilogit(psi.int[1,si]))

  N.br[ce,1,si] ~ dpois(exp(lamb.int[1,si])) # * z.br[ce,si])

  for (yr in 2:n.yr) {
    S.br[ce,yr,si] ~ dbin(prod(phi[yr-1,,si]) * (1 - ilogit(eps.int[1,si])), N.br[ce,yr-1,si])   

 ###  number of effective breeders:
    b.br[ce,yr,si] <- N.br[ce,yr-1,si] * (1 - ilogit(eps.int[1,si])) + (exp(iota.int[1,si]) / prod(phi[yr-1,,si]))
    
    G.br[ce,yr,si] ~ dpois(rho[yr-1,si] * prod(phi[yr-1,2:4,si]) * b.br[ce,yr,si])
    
    I.br[ce,yr,si] ~ dpois(exp(iota.int[1,si])) #dpois(z.br[ce,si] * exp(iota.int[1,si]))

    N.br[ce,yr,si] <- S.br[ce,yr,si] + G.br[ce,yr,si] + I.br[ce,yr,si]

  }  ## END yr loop

}  ##  END site loop


######################################################################################
####  Winter dynamics:  ##############################################################


for (tr in 1:n.sites[2]) {

#  z.nb[tr,si] ~ dbern(ilogit(psi.int[2,si]))
  
  N.nb[tr,1,si] ~ dpois(exp(lamb.int[2,si])) # * z.nb[tr,si])
  
  for (yr in 2:n.yr) {

    S.nb[tr,yr,si] ~ dbin(prod(phi[yr-1,3:4,si])*prod(phi[yr,1:2,si])*(1 - ilogit(eps.int[2,si])), N.nb[tr,yr-1,si])   ##  "survivors from year yr-1 to yr
 ###  effective number of breeders:
  b.nb[tr,yr,si] <- N.nb[tr,yr-1,si] * prod(phi[yr-1,3:4,si]) * (1 - ilogit(eps.int[2,si])) + (exp(iota.int[2,si]) / prod(phi[yr,1:2,si]))
  
    G.nb[tr,yr,si] ~ dpois(rho[yr,si]/2 * phi[yr,2,si] * b.nb[tr,yr,si])
    
    I.nb[tr,yr,si] ~ dpois(exp(iota.int[2,si]))  #dpois(z.nb[tr,si] * exp(iota.int[2,si]))

    N.nb[tr,yr,si] <- S.nb[tr,yr,si] + G.nb[tr,yr,si] + I.nb[tr,yr,si]

  }  ## END yr loop

}  ##  END site loop


########################################################################################
###  Observation:  ##############################################################

for (yr in 1:n.yr) {
  for (ce in 1:n.sites[1]) {
    obs.br[ce,yr,si] ~ dbin(p[1], N.br[ce,yr,si]) 
  }

  for (tr in 1:n.sites[2]) {
    obs.nb[tr,yr,si] ~ dbin(p[2], N.nb[tr,yr,si]) 
  }
  
  N.br.mn[yr,si] <- mean(N.br[,yr,si])
  N.nb.mn[yr,si] <- mean(N.nb[,yr,si])

}

## average population growth:

  grwth.br[1:(n.yr-1),si] <- N.br.mn[2:n.yr,si]/N.br.mn[1:(n.yr-1),si]
  grwth.nb[1:(n.yr-1),si] <- N.nb.mn[2:n.yr,si]/N.nb.mn[1:(n.yr-1),si]

  grwth.br.mn[si] <- prod(grwth.br[,si])^(1/(n.yr-1))   ## geometric mean
  grwth.nb.mn[si] <- prod(grwth.nb[,si])^(1/(n.yr-1))   ## geometric mean

} ## END loop over sims

}  ## END model
