##  Simulation master, FAC-IPM.
##
##
##  14 jul 2022; cleaned up 31 may 2024
#######################################

##  this should be called with Rscript.


##########################################
## General plan, for simulation exercise: 
## There are two entities involved, for each simulation exercise: a generator, and a JAGS model.

## 1. Call the generator, which gives the `dat` and `true` lists.  The `dat is fed to the JAGS model.
## 2. Fit the desired JAGS model.
## 3. Extract posteriors.
## 4. Repeat 1-3 many times.
## 5. Compare each iteration's posteriors with priors / true values.  Summarize the set of comparisons.
##########################################


##########################################
## Contents of this script:
##  Block 1: ingest and name Rscript arguments
##        2: simulate datasets
##        3: set up JAGS models
##        4: run JAGS models
##########################################



#######################################################################
##  BEGIN Block 1: ingest and name Rscript arguments                ##
#####################################################################

##  arguments:
##  1. "par.file" -- filename for R file holding the `par.lst` that gives the simulation's parameter values (no suffix -- should be '.r')
##  2. "mod.file" -- filename for JAGS model to be fit (no suffix -- should be '.r')
##  3. "ZIP" -- should there be zero inflation, in simulations? one of c("yesZ", "noZ")
##  4. "n.sim" -- number of simulations to put into this chunk
##  5. "br.n.site" -- number of breeding season sites (default 100)
##  6. "nb.n.site" -- number of winter sites (default 100)
##  7. "ann.var" -- Should there be annual variation in phi and rho?  one of c("yesAnn", "noAnn")
##  8. "chunk" -- chunk ID

args <- commandArgs(trailingOnly = TRUE)

names(args) <- c("par.file", "mod.file", "ZIP", "n.sim", "br.n.site", "nb.n.site", "ann.var", "chunk")

stopifnot(args["ZIP"] %in% c("noZ", "yesZ"))
stopifnot(args["ann.var"] %in% c("yesAnn", "noAnn"))

#####################################################################
##  END Block 1                                                     ##
#######################################################################


#######################################################################
##  BEGIN Block 2: simulate datasets                                ##
#####################################################################

##  Load parameters for the generating model:
load(paste0('sims/', args[1], '.r'))

##  Load the data-generating function and run it to create a number of datasets:
source('functions/effBreed_dec2023.r')

if (args["ann.var"] == "yesAnn") {
  sim.lst <- gen.noCovs(n.sim = as.integer(args["n.sim"]), n.sites = as.integer(args[c("br.n.site","nb.n.site")]), psi.mn = if(args["ZIP"] == "yesZ") c(0.7,0.7) else NULL)
} else {
  sim.lst <- gen.noCovs(n.sim = as.integer(args["n.sim"]), n.sites = as.integer(args[c("br.n.site","nb.n.site")]), psi.mn = if(args["ZIP"] == "yesZ") c(0.7,0.7) else NULL, phi.sd = c(0,0,0,0), rho.sd = 0)
}

##  Find the number of separate datasets created; save that in a list element:
sim.lst <- within(sim.lst, 
  dat$n.sim <- dim(dat$obs.br)[3]
)

#####################################################################
##  END Block 2                                                     ##
#######################################################################



#######################################################################
##  BEGIN Block 3: set up JAGS models                               ##
#####################################################################

##  Make initial values for latent, integer-valued state variables, for use in JAGS model run:
source('functions/init_maker_simDim.r')
init.lst <- init.maker(sim.lst$dat, simDim = T)

##  Identify variables to monitor:
vars <- c(
  "phi.int", "phi.sd",
  "rho.int", "rho.sd",
  "psi.int", 
  "eps.int", "eps.sd",
  "iota.int", "iota.sd",
  "lamb.int", 
  "grwth.br.mn", "grwth.nb.mn")


##  Load functions to run JAGS models in parallel:
source('../fit/functions/rjags_in_parallel.r')

library(parallel)
library(rjags)
load.module("lecuyer")

n.cores <- 6
rnd <- parallel.seeds("lecuyer::RngStream", n.cores)

dat <- sim.lst$dat

jl.unit <- list(iter.burn = 1e4,
                iter.samp = 1e4,
                var = vars,
                dat = dat, #c(dat, list(rem.p = inits.stuff$rem.v)),
                thin = 5,
                adapt = 1e3,
                file.nm = paste0('jags/', args["mod.file"], '.r'),
                inits = init.lst,    ### inits.stuff.proj
                drop = F)

jl <- lapply(rnd, FUN=c, jl.unit)                 

#####################################################################
##  END Block 3                                                     ##
#######################################################################


#######################################################################
##  BEGIN Block 4: run JAGS models                                  ##
#####################################################################

time <- system.time(
  out.pre <- mclapply(jl, FUN=jags.fun, mc.cores = n.cores)  ## mclapply() parallelizes this!
)

out <- mcmc.list(sapply(out.pre, FUN = function (l) { as.mcmc(get("js", l)) }, simplify=F))
mod <- sapply(out.pre, FUN = function (l) { get("jm", l) }, simplify=F)

fnm <- paste0(args["chunk"], "_", args["par.file"], "_", args["mod.file"], "_", Sys.Date())
ts <- Sys.time()

if (!dir.exists('output/mcmc')) {
  dir.create('output/mcmc')
}

save(out, mod, sim.lst, time, ts, file = paste0('output/mcmc/', fnm, '.R'))

#####################################################################
##  END Block 4                                                     ##
#######################################################################



