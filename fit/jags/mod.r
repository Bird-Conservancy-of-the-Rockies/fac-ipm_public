##  Parallel Dail-Madsen models, with demographic parameters informed by nest monitoring and telemetry data.
##
##  Supporting information for manuscript 'Linked dynamic abundance models as the basis for a full-annual-cycle, integrated population model'
##
##  04 dec 2024
###################################


data {
  br.B <- max(br.obs.x)
  br.n.obs <- length(br.obs.x)
  
  nb.B <- max(nb.obs.x)
  nb.n.obs <- length(nb.obs.x)

  br.n.ct <- length(br.ct.yr)
  nb.n.ct <- length(nb.ct.yr)

## Standard log-linear and logit-linear priors for normal distribution precision:  
  log.tau <- 0.15     ## altered from 0.01, 05 jul 2023
  logit.tau <- 0.37   ## Jeffreys prior
  
## For use with half-t priors: k, the degrees of freedom.  Settled on these by messing around with values, after seeing Plummer's note that k = 4 is a common value, in the JAGS manual.
  log.k <- 10
  logit.k <- 5
  
#  log.lapl <- 1/2.582         ## for Laplace priors (i.e., ddexp()), where ddexp "tau" = 1/sigma, where sigma is kind of like SD of normal distr.
#  logit.lapl <- 1/1.644       ## for Laplace priors (i.e., ddexp()), where ddexp "tau" = 1/sigma, where sigma is kind of like SD of normal distr.
  
#  small <- 0.001

## when using the 'weak philopatry' dynamics in breeding season, this needs to be '1': otherwise, '4'
  first.br.yr <- 5
  first.nb.yr <- 1

## mean of all counts:
  mn.br.ct <- mean(br.ct.tot)
  mn.nb.ct <- mean(nb.ct.tot)

## zeros for annual random effects on phi (four adult survivals, plus offspring during breeding season)
  phi.mu.0 <- c(0,0,0,0,0)
## ...and movement:  
  ie.mu.0 <- c(0,0,0,0)  

## dimensions:
  dbr.y.off.dim <- dim(dbr.y.off)
  dbr.y.ad.dim <- dim(dbr.y.ad)

  nstPer <- 21
  fldgPer <- 14
## in this version of the model, trap period may be defined here:
  trapPer <- 8

## predictor coefficients:
##  Check to make sure these are right, if you see, "Unable to resolve relations"!
#  n.phi.br.covs <- length(phi.br.covs)
  n.phi.off.covs <- length(phi.off.covs)
#  n.psi.br.covs <- length(psi.br.covs)   
  n.lamb.br.covs <- length(lamb.br.covs)
  n.eps.br.covs <- length(eps.br.covs)
  n.iota.br.covs <- length(iota.br.covs)
  
  n.phi.nb.covs <- length(phi.nb.covs)
#  n.psi.nb.covs <- length(psi.nb.covs)
  n.lamb.nb.covs <- length(lamb.nb.covs)
  n.eps.nb.covs <- length(eps.nb.covs)
  n.iota.nb.covs <- length(iota.nb.covs)

## total dimensions of covariate arrays:
  dim.br.x <- dim(br.x)
  dim.nb.x <- dim(nb.x)

## number of SRV plots, years:
  n.srv.plot <- length(nb.plot.n.tran)
  n.srv.yr <- max(dnb.yr)

}  ## END data


model {

### Annual survival parameters:
for (i in 1:5) {  ## the four adult survivals, plus offspring, #5.
#  phi.int[i] ~ ddexp(0, logit.tau*2)  ## times two, because the ddexp() variance is twice that of a normal
## SD for annual random effect: hyperparameter should help with forecasting.
#  phi.sd[i] ~ dt(0,1,logit.k) T(0,) #dunif(0, 1.64399)  ## 1.64399 = logit.tau^(-0.5)
  phi.sd[i] ~ dunif(0, 0.75)
}

phi.int[1] ~ ddexp(0, logit.tau*2)  ## times two, because the ddexp() variance is twice that of a normal
phi.mn.2 ~ dbeta(3,3)
phi.int[2] <- logit(phi.mn.2)
phi.int[3] ~ ddexp(0, logit.tau*2)  ## times two, because the ddexp() variance is twice that of a normal
phi.mn.4 ~ dbeta(3,3)
phi.int[4] <- logit(phi.mn.4)
phi.int[5] ~ ddexp(0, logit.tau*2)  ## times two, because the ddexp() variance is twice that of a normal


############################################################################
##  Linear predictor intercepts: it's just easier to specify these by name here.

eps.br.int ~ dnorm(0,logit.tau)
ie.sd[1] ~ dunif(0, 0.75)   ## breeding emigration random effect SD

lamb.br.int ~ dnorm(0,log.tau)
psi.br.int ~ dnorm(0,logit.tau)
iota.br.int ~ dnorm(0,log.tau)
ie.sd[2] ~ dunif(0, 0.75)   ## breeding immigration random effect SD

eps.nb.int ~ dnorm(0,logit.tau)
ie.sd[3] ~ dunif(0, 0.75)   ## winter emigration random effect SD

lamb.nb.int ~ dnorm(0,log.tau)
psi.nb.int ~ dnorm(0,logit.tau)

iota.nb.int ~ dnorm(0,log.tau)
ie.sd[4] ~ dunif(0, 0.75)   ## winter immigration random effect SD


##########   offspring season-wise survival, for count model:
beta.nst ~ dnorm(0,logit.tau)
beta.fldg ~ dnorm(0,logit.tau)
beta.juv <- 0   ## to make the 

##########   adult trap effect on survivals: these are effects on DAILY survival!
beta.br.trap ~ dnorm(0,logit.tau)
beta.nb.trap ~ dnorm(0,logit.tau)

##########   suppose there's some systemic discrepancy between adult survival, as estimated from the count data, versus from the telemetry data: these are effects on SEASONAL survivals!
cmr.br <- 0  #~ dnorm(0,logit.tau)
cmr.nb <- 0  #~ dnorm(0,logit.tau)

##########   breeding season adult survival: no covariates in this model version
#for (i in 1:n.phi.br.covs) {
#  beta.phi.br[i] ~ dnorm(0, logit.tau)
#}

##########   breeding season offspring survival:
for (i in 1:n.phi.off.covs) {
  beta.phi.off[i] ~ dnorm(0, logit.tau)
}


##########   breeding season initial abundance:
for (i in 1:n.lamb.br.covs) {
  beta.lamb.br[i] ~ dnorm(0, log.tau)
}


##########   breeding season emigration:
for (i in 1:n.eps.br.covs) {
  beta.eps.br[i] ~ dnorm(0, logit.tau)
}


##########   breeding season immigration:
for (i in 1:n.iota.br.covs) {
  beta.iota.br[i] ~ dnorm(0, log.tau)
}


##########   winter initial abundance:
for (i in 1:n.lamb.nb.covs) {
  beta.lamb.nb[i] ~ dnorm(0, log.tau)
}


##########   winter adult survival:
for (i in 1:n.phi.nb.covs) {
  beta.phi.nb[i] ~ dnorm(0, logit.tau)
}


##########   winter emigration:
for (i in 1:n.eps.nb.covs) {
  beta.eps.nb[i] ~ dnorm(0, logit.tau)
}


##########   winter immigration:
for (i in 1:n.iota.nb.covs) {
  beta.iota.nb[i] ~ dnorm(0, log.tau)
}


##  How long is each season?
##  Must know this in order to convert seasonal survival to daily survival...
##  For now, assume they are all equal.
seas.len[1] <- 365/4        ## summer
seas.len[2] <- seas.len[1]  ## fall
seas.len[3] <- seas.len[1]  ## winter
seas.len[4] <- seas.len[1]  ## spring

##  Average length of the juvenile phase is the length of summer, minus the average day-of-summer when nests fledged.
juvPer <- seas.len[1] - dbr.mn.fldg.day


####  Year-level random effects on survival rates:
## I could estimate correlation between annual effects on survival rates:
##   the four seasonal adult survivals, and offspring survival during breeding season.
##   should I do all potential correlations?  (5^2 - 5)/2 = 10...

for (i in 1:5) {
  phi.Sigma[i,i] <- phi.sd[i]^2
}

for (i in 2:5) {
  phi.cor[i,1] <- 0 #~ dunif(-1,1)
  phi.Sigma[1,i] <- phi.cor[i,1] * phi.sd[1] * phi.sd[i]
  phi.Sigma[i,1] <- phi.Sigma[1,i]
}

for (i in 3:5) {
  phi.cor[i,2] <- 0 #~ dunif(-1,1)
  phi.Sigma[2,i] <- phi.cor[i,2] * phi.sd[2] * phi.sd[i]
  phi.Sigma[i,2] <- phi.Sigma[2,i]
}

for (i in 4:5) {
  phi.cor[i,3] <- 0 #~ dunif(-1,1)
  phi.Sigma[3,i] <- phi.cor[i,3] * phi.sd[3] * phi.sd[i]
  phi.Sigma[i,3] <- phi.Sigma[3,i]
}

phi.cor[5,4] <- 0 #~ dunif(-1,1)
phi.Sigma[4,5] <- phi.cor[5,4] * phi.sd[4] * phi.sd[5]
phi.Sigma[5,4] <- phi.Sigma[4,5]



####  Year-level random effects on movement rates:  this is basically allowing movement out of the monitored the population; or into it from without.
##  I could estimate correlation between annual effects on survival rates:
##   the four seasonal adult survivals, and offspring survival during breeding season.
##   should I do all potential correlations?  (4^2 - 4)/2 = 6...


for (i in 1:4) {
  ie.Sigma[i,i] <- ie.sd[i]^2
}

for (i in 2:4) {
  ie.cor[i,1] <- 0 #~ dunif(-1,1)
  ie.Sigma[1,i] <- ie.cor[i,1] * ie.sd[1] * ie.sd[i]
  ie.Sigma[i,1] <- ie.Sigma[1,i]
}

for (i in 3:4) {
  ie.cor[i,2] <- 0 #~ dunif(-1,1)
  ie.Sigma[2,i] <- ie.cor[i,2] * ie.sd[2] * ie.sd[i]
  ie.Sigma[i,2] <- ie.Sigma[2,i]
}

ie.cor[4,3] <- 0 #~ dunif(-1,1)
ie.Sigma[3,4] <- ie.cor[4,3] * ie.sd[3] * ie.sd[4]
ie.Sigma[4,3] <- ie.Sigma[3,4]

cltch.sd ~ dt(0,1,logit.k) T(0,)   ## Half-Cauchy distribution.
###  Random effect of clutch on offspring survival:
for (i in 1:dbr.n.nst) {
  alpha.cltch[i] ~ dnorm(0, cltch.sd^(-2))
}


###  Clutch size is universal: beta-binomial?
cl.b1 ~ dgamma(1,1)
cl.b2 ~ dgamma(1,1)
cl.pr ~ dbeta(cl.b1, cl.b2)

for (n in 1:dbr.n.cl) {
  dbr.cl[n] ~ dbinom(cl.pr, 6)
}

##  Expected clutch size.
cl.expect <- cl.pr * 6

###  Breeding demographic sites: 'NGP'__________________________________________________________
for (o in 1:dbr.y.off.dim[1]) {    ## BEGIN `o` loop through offspring

###  FUNCTION : Offspring daily survival at demographic study sites

## These are just the definitions from below in the 'cells' loop:  
## phi.d.off is the "intercept" daily survival rate: this means that the overall mean daily survival rate MIGHT DIFFER from the mean rate according to the count model.
##  to ensure that the mean rates are equal, I guess phi.d.juv's prior could be replaced by a sum-to-zero contraint: phi.d.juv = -(phi.d.nst + phi.d.fldg)
## convert to daily survival: phi.d.off---since it is spread over the entire breeding season and not all nests are initiated at the beginning of breeding season---might be interpreted as the survival probability of a "theoretical egg" that may or may not be laid: such a "hypothetical ovum" must survive until conception, through nest building, laying, etc.  So it doesn't matter when in the season the clutch was actually laid.
  phi.d.off.int[o] <- logit( phi.off.ce[br.plot2cell[dbr.off.pl[o]], yr.dbr2un[dbr.off.yr[o]]]^(1/seas.len[1]) ) + alpha.cltch[dbr.off.nst[o]]
  logit(phi.d.nst[o]) <- phi.d.off.int[o] + beta.nst
  logit(phi.d.fldg[o]) <- phi.d.off.int[o] + beta.fldg
  logit(phi.d.juv[o]) <- phi.d.off.int[o]   ## juvenile survival is the 'reference condition' for offspring survival

## for checking:
#  phi.off.diff[o] <- phi.off.ce[br.plot2cell[dbr.off.pl[o]], yr.dbr2un[dbr.off.yr[o]]] - (phi.nst[o] * phi.fldg[o] * phi.juv[o]) 
  
  for (d in 1:(dbr.k.off[o]-1)) {  ## post-trap Period
###  FUNCTION : Assign survival values to days (of life-after-capture)
    phi.d.off[o,d] <- ifelse( d <= nstPer,
                             phi.d.nst[o],
                             ifelse(d <= (nstPer + fldgPer),
                               phi.d.fldg[o],
                               phi.d.juv[o]) )
####  LIKELIHOOD : Offspring daily survival at demographic study sites                             
    dbr.y.off[o,d+1] ~ dbern( phi.d.off[o,d] * dbr.y.off[o,d] )
####  For posterior predictive check:
#    dbr.y.off.new[o,d+1] ~ dbern( phi.d.off[o,d] * dbr.y.off[o,d] )    
#    E.dbr.y.off[o,d+1] <- phi.d.off[o,d] * dbr.y.off[o,d]
    
  }

}    ## END `o` loop through offspring



for (a in 1:dbr.y.ad.dim[1]) {     ## BEGIN `a` loop through telemetered breeding adults

###  FUNCTION : Breeding adult (seasonal) survival at demographic study sites
###    IN this 'parShare' version, just take the survival values from the appropriate 'grid cells' that represent the demographic sites:
### the base (average) daily survival:
  phi.d.dbr.mn[a] <- ilogit(logit(phi.br[br.plot2cell[dbr.ad.pl[a]], yr.dbr2un[dbr.ad.yr[a]]]) + cmr.br)^(1/(seas.len[1]))

  for (d in 1:(dbr.k.ad[a]-1)) {  ## post-trap Period
###  FUNCTION : Assign survival values to days (of life-after-capture)
    phi.d.dbr[a,d] <- ifelse(d <= trapPer,
                             ilogit(logit(phi.d.dbr.mn[a]) + beta.br.trap),
                             phi.d.dbr.mn[a])
####  LIKELIHOOD : Breeding adult daily survival at demographic study sites                             
    dbr.y.ad[a,d+1] ~ dbern( phi.d.dbr[a,d] * dbr.y.ad[a,d] )
  }

}   ## END `a` loop over individuals



###  Winter telemetry sites: 'SRV'__________________________________________________________
### First order of business is to summarize survival rate across transects within each plot:

for (yr in 1:n.srv.yr) {
  for (pl in 1:n.srv.plot) {
    for (tr in 1:nb.plot.n.tran[pl]) {
      phi.dnb.tr[yr,pl,tr] <- phi.nb[nb.plot2tran[pl,tr], yr.dnb2un[yr]]   ## conversion of plot and year indexing, from count to telemetry datasets.
    }
    phi.dnb.pl[yr,pl] <- mean(phi.dnb.tr[yr, pl, 1:nb.plot.n.tran[pl]])
  }
}

for (a in 1:dnb.n.ind) {     ## BEGIN `a` loop through telemetered winter adults

###  FUNCTION : Winter adult daily survival at demographic study sites
### the base (average) daily survival:
  phi.d.dnb.mn[a] <- ilogit(logit(phi.dnb.pl[dnb.yr[a], dnb.site[a]]) + cmr.nb)^(1/(seas.len[3]))
  
####  LIKELIHOOD : Winter adult daily survival at demographic study sites
##  Condensing things here:
  for (d in 1:(dnb.k[a]-1)) {  ## post-trap Period
###  FUNCTION : Assign survival values to days (of life-after-capture)
    phi.d.dnb[a,d] <- ifelse( d <= trapPer,
                             ilogit(logit(phi.d.dnb.mn[a]) + beta.nb.trap),
                             phi.d.dnb.mn[a] )
                             
    dnb.y[a,d+1] ~ dbern( phi.d.dnb[a,d] * dnb.y[a,d] )
  }

}     ## END `a` loop through telemetered winter adults





### Productivity and survival rates at sites: this is per male(ish)!
###  Fit to demography sites here, first:
for (yr in 1:n.yr) {

  alpha.phi[yr,1:5] ~ dmnorm.vcov(phi.mu.0, phi.Sigma)
  alpha.ei[yr,1:4] ~ dmnorm.vcov(ie.mu.0, ie.Sigma)

## Migratory season survivals:
  ## Autumn survival:
  logit(phi[yr,2]) <- phi.int[2] + alpha.phi[yr,2]

  ## Spring daily survival:s
  logit(phi[yr,4]) <- phi.int[4] + alpha.phi[yr,4]
  
  
## Breeding:
  for (ce in 1:br.n.cell) {    ## BEGIN `ce` loop over grid cells
  
## Don't need phi.off for every cell: only for demography site-yr combinations!    
##  Instead, need to calculate nst.phi and juv.phi:

###  FUNCTION : offspring daily survival, at grid cells
#    logit(phi.off.ce[ce,yr]) <- phi.int[5] +           
#                                  alpha.phi[yr,5] +                           ## year random effect
#                                  beta.phi.off %*% br.x[ce, yr, phi.off.covs]  ## covariate fixed effects

###  FUNCTION : nest daily survival, at grid cells
###  Next 5 statements could be condensed into one!
#    logit(phi.nst.ce[ce,yr]) <- phi.int[5] +           
#                                  alpha.phi[yr,5] +                           ## year random effect
#                                  beta.phi.nst %*% br.x[ce, yr, phi.off.covs]  ## covariate fixed effects

#    logit(phi.fldg.ce[ce,yr]) <- logit(phi.nst.ce[ce,yr]) + beta.phi.fldg

#    logit(phi.juv.ce[ce,yr]) <- logit(phi.nst.ce[ce,yr]) + beta.phi.juv

#    phi.off.ce[ce,yr] <- phi.nst.ce[ce,yr] * phi.fldg.ce[ce,yr] * phi.juv.ce[ce,yr]
#   

#    rho.ce[ce,yr] <- (phi.off.ce[ce,yr] * cl.expect)     

####  FUNCTIONs : seasonal productivity, at grid cells.
##  These are "period-wise" survivals (not daily).
logit(phi.off.ce[ce,yr]) <- phi.int[5] + beta.phi.off %*% br.x[ce, yr, phi.off.covs] + alpha.phi[yr,5]
#phi.fldg.ce[ce,yr] <- ilogit(logit(phi.nst.ce[ce,yr]^(1/nstPer)) + beta.phi.fldg)^(fldgPer)     ## beta.phi.fldg is now the effect on DAILY survival of being a fledgling, vs. an egg/nestling
#phi.juv.ce[ce,yr] <- ilogit(logit(phi.nst.ce[ce,yr]^(1/nstPer)) + beta.phi.juv)^(seas.len[1] - (nstPer + fldgPer))   ## beta.phi.juv is now the effect on DAILY survival of being a juvenile, vs. an egg/nestling

####  FUNCTIONs : seasonal productivity, at grid cells
rho.ce[ce,yr] <- cl.expect * phi.off.ce[ce,yr]     ## essentially per male reproduction!


###  FUNCTIONs : breeding adult survival, at grid cells
    logit(phi.br[ce,yr]) <- phi.int[1] + 
                               alpha.phi[yr,1] #+              ## year random effect
#                               beta.phi.br %*% br.x[ce,yr,phi.br.covs]
  
  }   ## END `ce` loop over cells


###  Winter:
  for (tr in 1:nb.n.tran) {    ## BEGIN `tr` loop over transects
  
###  FUNCTION : winter adult survival, at transects
    logit(phi.nb[tr,yr]) <- phi.int[3] + 
                               alpha.phi[yr,3] +                  ## year random effect
                               beta.phi.nb %*% nb.x[tr,yr,phi.nb.covs]
  
  }  ## END `tr` loop over transects


###  FUNCTIONs : Average seasonal adult survivals, and productivity:
  logit(phi[yr,1]) <- phi.int[1] #+ beta.phi.br %*% br.x.mn.calc[yr,phi.br.covs] + alpha.phi[yr,1] ## breeding adult
  logit(phi[yr,3]) <- phi.int[3] + beta.phi.nb %*% nb.x.mn.calc[yr,phi.nb.covs] + alpha.phi[yr,3] ## winter adult
  logit(phi[yr,5]) <- phi.int[5] + beta.phi.off %*% br.x.mn.calc[yr, phi.off.covs] + alpha.phi[yr,5]   ## offspring

  logit(eps.br.mn[yr]) <- eps.br.int + beta.eps.br %*% br.x.mn.calc[yr, eps.br.covs] + alpha.ei[yr,1]   ## net breeding epsilon
  logit(iota.br.mn[yr]) <- iota.br.int + beta.iota.br %*% br.x.mn.calc[yr, iota.br.covs] + alpha.ei[yr,2]   ## new breeding iota

  logit(eps.nb.mn[yr]) <- eps.nb.int + beta.eps.nb %*% nb.x.mn.calc[yr, eps.nb.covs] + alpha.ei[yr,3]   ## net breeding epsilon
  logit(iota.nb.mn[yr]) <- iota.nb.int + beta.iota.nb %*% nb.x.mn.calc[yr, iota.nb.covs] + alpha.ei[yr,4]   ## new breeding iota

  rho[yr] <- phi[yr,5] * cl.expect     ## this is per male reproduction, essentially!

###  FUNCTIONs : Average annual survival (pre-breeding census: beginning of summer) 
  phi.ann[yr] <- prod(phi[yr,1:4])

###  FUNCTIONs : Average annual growth (pre-breeding census: beginning of summer) 
## this is for yr --> yr+1  !!
  grwth.pop[yr] <- phi.ann[yr] +
                   (prod(phi[yr,2:4]) * rho[yr]/2)  ## the '2' translates per male to per capita reproduction!

}  ## END `yr` loop over years


########################################
##  Psi coefficients never subject to selection:

for (i in 1:3) {
  beta.psi.br[i] ~ dnorm(0, logit.tau)   ## note the "times two"!
}

for (ce in 1:br.n.cell) {
  logit(psi.br[ce]) <- psi.br.int +
                       beta.psi.br[1] * br.coord[ce,1] +                 ## longitude
                       beta.psi.br[2] * br.coord[ce,2] +                 ## latitude
                       beta.psi.br[3] * br.coord[ce,1] * br.coord[ce,2]  ## interaction
}


for (yr in 1:n.yr) {
###  Maybe the "mean covariate value" that's used to calculate population-wide survivals should be the mean across *occupied* sites: "psi.br %*% br.x" and  "psi.nb %*% nb.x"
  br.x.mn.calc[yr,1:dim.br.x[3]] <- (psi.br/sum(psi.br)) %*% br.x[,yr,] 
  nb.x.mn.calc[yr,1:dim.nb.x[3]] <- (psi.nb/sum(psi.nb)) %*% nb.x[,yr,]
}


####  RANDOM EFFECTS of NB stratum: non-breeding rates  #######
psi.nb.sd ~ dunif(0,5)  ## logit

for (g in 1:nb.n.strat) {
  alpha.psi.nb[g] ~ dnorm(0, psi.nb.sd^(-2))
}

########################################################


########################################
##  Psi NB:

for (tr in 1:nb.n.tran) {
  logit(psi.nb[tr]) <- psi.nb.int + alpha.psi.nb[nb.tr.strat[tr]]
}


for (yr in 2:n.yr) {

#########   Population-level recruitment: these annual intercepts will be put on the log scale when they are used, downstream.
  gamma.br.mn[yr] <- rho[yr-1] * prod(phi[yr-1,2:4])
  gamma.nb.mn[yr] <- rho[yr]/2 * phi[yr,2]       ### N.B. - the `2` here is to translate per male reproduction to per capita! 

#########   Population-level survival:
  phiA.br.mn[yr] <- prod(phi[yr-1,1:4])
  phiA.nb.mn[yr] <- prod(phi[yr-1,3:4]) * prod(phi[yr,1:2])
  
}


######################################################################################
####  Summer dynamics:  ##############################################################

for (ce in 1:br.n.cell) {
  
  log(lambda.br[ce]) <- lamb.br.int + beta.lamb.br %*% br.x[ce,first.br.yr,lamb.br.covs]

  z.br[ce] ~ dbern(psi.br[ce])
  N.br[ce,first.br.yr] ~ dpois(lambda.br[ce] * z.br[ce])
  
  for (yr in (first.br.yr+1):n.yr) {
###  Here's where the site meets the FAC!
    phiA.br[ce,yr] <- phi.br[ce,yr-1] * prod(phi[yr-1,2:4])

##  Here's where the site meets the FAC!
    logit(epsilon.br[ce,yr]) <- eps.br.int + alpha.ei[yr,1] + beta.eps.br %*% br.x[ce,yr,eps.br.covs]

    gamma.br[ce,yr] <- rho.ce[ce,yr-1] * prod(phi[yr-1,2:4])
    
    log(iota.br[ce,yr]) <- iota.br.int +
                           alpha.ei[yr,2] + beta.iota.br %*% br.x[ce,yr,iota.br.covs]

    S.br[ce,yr] ~ dbin(phiA.br[ce,yr] * (1 - epsilon.br[ce,yr]), N.br[ce,yr-1])   ##  "survivors from year yr-1 to yr
    b.br[ce,yr] <- N.br[ce,yr-1] * (1 - epsilon.br[ce,yr]) + (z.br[ce] * iota.br[ce,yr])/ prod(phi[yr-1,])  ## effective number of breeders contributing offspring to this site
    G.br[ce,yr] ~ dpois(gamma.br[ce,yr] * b.br[ce,yr])
    I.br[ce,yr] ~ dpois(z.br[ce] * iota.br[ce,yr])  
    N.br[ce,yr] <- S.br[ce,yr] + G.br[ce,yr] + I.br[ce,yr]

  }  ## END yr loop
}  ##  END site loop


######################################################################################
####  Winter dynamics:  ##############################################################

for (tr in 1:nb.n.tran) {

  log(lambda.nb[tr]) <- lamb.nb.int + beta.lamb.nb %*% nb.x[tr,first.nb.yr,lamb.nb.covs]
  z.nb[tr] ~ dbern(psi.nb[tr])
  N.nb[tr,first.nb.yr] ~ dpois(lambda.nb[tr] * z.nb[tr])

  for (yr in (first.nb.yr+1):n.yr) {

###  There is no separate 'recruitment' process: juveniles are 'attached' to adults in the fall; emigration occurs as all birds, SY and ASY, show up on the wintering grounds.

###  Here's where the site meets the FAC!
    phiA.nb[tr,yr] <- phi.nb[tr,yr-1] * phi[yr-1,4] * prod(phi[yr,1:2])

    logit(epsilon.nb[tr,yr]) <- eps.nb.int + alpha.ei[yr,3] + beta.eps.nb %*% nb.x[tr,yr,eps.nb.covs]
     
    log(iota.nb[tr,yr]) <- iota.nb.int + alpha.ei[yr,4] + beta.iota.nb %*% nb.x[tr,yr,iota.nb.covs]
    
    S.nb[tr,yr] ~ dbin(phiA.nb[tr,yr] * (1 - epsilon.nb[tr,yr]), N.nb[tr,yr-1])                ##  "survivors from year yr-1 to year yr"
    b.nb[tr,yr] <- N.nb[tr,yr-1] * prod(phi[yr-1,3:4]) * (1 - epsilon.nb[tr,yr]) + (z.nb[tr] * iota.nb[tr,yr])/prod(phi[yr,1:2])  ## effective number of breeders contributing recruits to this site
    G.nb[tr,yr] ~ dpois(gamma.nb.mn[yr] * b.nb[tr,yr])
    I.nb[tr,yr] ~ dpois(z.nb[tr] * iota.nb[tr,yr])
    N.nb[tr,yr] <- S.nb[tr,yr] + G.nb[tr,yr] + I.nb[tr,yr]

  }  ##  END yr loop
}  ##  END site loop


#########################################################################################
####  Summer observation:  ##############################################################

br.sigma.int ~ dnorm(0, log.tau)
br.sigma.sd ~ dt(0, 1, log.k) T(0,)

for (yr in 4:n.yr) {    ##  there are no IMBCR surveys for years 1,2,3
  br.sigma.alpha[yr] ~ dnorm(0, br.sigma.sd^(-2))
  log(br.sigma[yr]) <- br.sigma.int + br.sigma.alpha[yr]

  br.denom[yr] <- -(br.sigma[yr]^2) * (exp( -(br.B^2)/(2*br.sigma[yr]^2) ) - 1)
  nu[yr] <- 2 * pi * br.denom[yr]
  br.p[yr] <- nu[yr] / (pi*br.B^2)

  for (db in 1:br.n.bin) {
  
## pi.c works out to: ( g(c_j-1) - g(c_j) ) / (1 - g(B))
    br.pi.c[db,yr] <- ( exp(-(br.cut[db]^2)/(2*br.sigma[yr]^2)) - exp(-(br.cut[db+1]^2)/(2*br.sigma[yr]^2)) ) / (1 - exp(-(br.B^2)/(2*br.sigma[yr]^2)) )
  }  ## END db loop

  br.bin.tot[,yr] ~ dmulti(br.pi.c[,yr], br.yr.tot[yr])
}  ## END yr loop


for (k in 1:br.n.ct) {
##  Here, N.br is interpreted as the density in the full set of 16 point counts in an IMBCR grid cell.  If fewer than 16 counts were performed in a particular grid cell visit, br.ct.intens will be less than 1.0, and fewer individuals will potentially be available to be counted in the visited point counts.
  br.ct.avail[k] ~ dpois(N.br[br.ct.cell[k], br.ct.yr[k]] * br.ct.intens[k])    ## availability
  br.ct.tot[k] ~ dbin(br.p[br.ct.yr[k]], br.ct.avail[k])                        ## detection     
}   ## END k loop

#-------------------



#########################################################################################
####  Winter observation:  ##############################################################

nb.sigma.int ~ dnorm(0, log.tau)
nb.sigma.sd ~ dt(0, 1, log.k) T(0,)   ## Half-Cauchy distribution

for (yr in 1:(n.yr-2)) {   ##  there are no winter surveys for years 15 and 16
  nb.sigma.alpha[yr] ~ dnorm(0, nb.sigma.sd^(-2))
  log(nb.sigma[yr]) <- nb.sigma.int + nb.sigma.alpha[yr]

  for (db in 1:nb.n.bin) {
        nb.p.db[db,yr] <- exp(-(nb.midpt[db]^2)/(2*nb.sigma[yr]^2)) * 1/nb.n.bin         ## Half-normal; "1/nb.n.bin" because all the distance bands have equal area.
        nb.pi.c[db,yr] <- nb.p.db[db,yr]/nb.p[yr]  #pi[k]/p.cap
  }  ## END db loop

  nb.p[yr] <- sum(nb.p.db[,yr]) #sum(pi[])

  nb.bin.tot[,yr] ~ dmulti(nb.pi.c[,yr], nb.yr.tot[yr])

}  ## END yr loop



for (k in 1:nb.n.ct) {
### The following adjustment by `rel.trans.len` accounts for the fact that transects have different lengths.  PLMX transects are only 500m, whereas GPCA transects are 1 km.  The interpretation of a PLMX 'site', then, is a 1-km transect, only half of which was surveyed!
  nb.ct.avail[k] ~ dpois(N.nb[nb.ct.tran[k],nb.ct.yr[k]] * nb.tr.km[k])     ## availability
  nb.ct.tot[k] ~ dbin(nb.p[nb.ct.yr[k]], nb.ct.avail[k])                    ## detection
}  ## END k loop


}  ## END model
