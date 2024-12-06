##  Master script to run JAGS model, with chains in parallel. 
##
##
##  05 dec 2024
############################################



##  This script is intended to be run from the commandline using Rscript.
##  For example, with working directory (pwd) set to 'fac-ipm_public/fit':
# Rscript scripts/run_mod.r 3 1e3 5e3 5e3 1 test

library(parallel)
library(rjags)
load.module("lecuyer")

source('functions/rjags_in_parallel.r')

args <- as.list(commandArgs(trailingOnly = TRUE))

names(args)<- c(
  "n.cores",   ##  JAGS parameters: Number of MCMC chains (to run in parallel), i.e., the number of cores.
  "n.adapt",   ##  JAGS parameters: Number of adaptation iterations.
  "n.burn",    ##  JAGS parameters: Number of burn-in iterations.
  "n.samp",    ##  JAGS parameters: Number of monitoring iterations.
  "n.thin",    ##  JAGS parameters: Thinning rate (inverse of it, actually).
  "label")     ##  extra label to make sure output filenames are unique

## if running interactively: example values
args <- list(
  n.cores = "3",
  n.adapt = "1e3",
  n.burn = "5e3",
  n.samp = "5e3",
  n.thin = "1",
  label = "test")



## convert some values to numbers, in fact most of them:
txt <- "label"
invisible(
  sapply(setdiff(names(args), txt), FUN = function (nm) {
    args[[nm]] <<- as.numeric(args[[nm]])
  })
)



# bring elements of args list into global environment, as named objects:
list2env(args, envir = .GlobalEnv)




##  bring in data object:
load('input/dat.r')

##  variables to monitor:
vars <- drop(read.csv("input/vars.csv", comment.char = "#", header = F)[,1])

##  bring in initial values:
load("input/inits.R")

##  set random number generators for cores: 
rnd <- parallel.seeds("lecuyer::RngStream", n.cores)


jl.unit <- list(adapt = n.adapt,
                iter.burn = n.burn,
                iter.samp = n.samp,
                var = vars,
                dat = dat, #c(dat, list(rem.p = inits.stuff$rem.v)),
                thin = n.thin,
                file.nm = paste0('jags/mod.r'),
                inits = init.lst,    ### inits.stuff.proj
                drop = F)

jl <- lapply(rnd, FUN=c, jl.unit)


t1 <- Sys.time()

####  run the model:
out.pre <- mclapply(jl, FUN=jags.fun, mc.cores = n.cores)

t2 <- Sys.time()
elap <- t2 - t1


####  separate MCMC output and model objects:
out <- mcmc.list(sapply(out.pre, FUN = function (l) { as.mcmc(get("js", l)) }, simplify=F))
mod <- sapply(out.pre, FUN = function (l) { get("jm", l) }, simplify=F)

## safety measure, in case subsequent steps to save output fail:
dir.create('jags_backup', recursive = TRUE)
save(out, mod, jl.unit, sc.lst, elap, file = paste0('~/jags_backup/jags_output_', Sys.Date(), '_', sample(x = 1e4, size = 1), '.R'))

###  this doesn't seem to work when there is a data block at the beginning of the JAGS model file: it won't re-run that portion!  Annoying...
#out.cont.pre <-  mclapply(mod, FUN=jags.cont, dat = dat, vars = vars, n.update = 1e4, n.monitor = 5e3, mc.cores = n.cores) 

#plot(out, ask = T)

  s <- summary(out, quantile = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975))

  compy <- system2("hostname", args = "-s", stdout = TRUE)
  fnm <- paste0(mod.nm, "_", Sys.Date(), "_", compy)

#save(out, mod, jl.unit, sc.lst, time, file = paste0('output/mcmc/', fnm, '.R'))

##  save a copy of the actual JAGS model code, along with everything else, since the only other way to do that would be to match the date of the run, which might not be the same as when the run was started:
  mod.code <- paste0(scan(jl.unit$file.nm, what = "character", sep = "\n"), collapse = "\n")

##  before attempting to save output, make sure the target directory exists:
  mcmc.path <- paste0('output/mcmc/', spp, '/ipm/')
  summ.path <- paste0('output/summ/', spp, '/ipm/')
  dir.create(mcmc.path, recursive = TRUE)
  dir.create(summ.path, recursive = TRUE)
  
  save(out, mod, jl.unit, sc.lst, elap, mod.code, file = paste0(mcmc.path, fnm, '_', label, '.R'))
  save(s, mod, jl.unit, sc.lst, elap, mod.code, file = paste0(summ.path, fnm, '_', label, '.R'))

