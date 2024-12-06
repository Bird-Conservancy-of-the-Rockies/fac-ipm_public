##  Collect results of fitting model to a number of simulated datasets.
##
##
##  26 jul 2022
############################################

gathrFits <- function (
  key,       ## vector of keywords to look for in files
  out.dir = 'output/mcmc',
  max.n.files = 10
) {

require(coda)
require(parallel)
require(abind)

source('functions/obj_maker_dimAdd.r')



lf <- list.files(out.dir)
tst.m <- sapply(key, FUN = function (s) { grepl(lf, pattern = s, fixed = F) })
tst.v <- apply(tst.m, MAR = 1, FUN = all)

fl <- lf[tst.v]
chunk <- as.integer(gsub(fl, pattern = "^([0-9]{1,2})\\_.*", replace = "\\1"))

fl <- fl[order(chunk)]

stopifnot(length(fl) == 10)

################################################   BLOCK 1 -- begin  ###########
#################    Retrieve true values
####
#   Let's see, all the 'truth' values should be the same, so just retrieve one such set:

load(paste0(out.dir, "/", fl[1]))

true <- sim.lst$true[!names(sim.lst$true) %in% c("psi", "rho", "phi")]

logit.par <- c("psi.int", "phi.int", "eps.int")
log.par <- c("lambda.int", "iota.int", "rho.int")

true.trans <- c(lapply(true[logit.par], FUN = plogis),
                lapply(true[log.par], FUN = exp))


#
####
#################
################################################   BLOCK 1 -- end  #############


################################################   BLOCK 2 -- begin  ###########
#################    Collect all MCMC results
####
#  In an array?  I guess I want to know the bias / coverage for each dataset...

want <- names(true)

fit.lst <- mclapply(1:length(fl), FUN = function (n) {
  load(paste0(out.dir, "/", fl[n]))
  sub.out <- as.mcmc.list(lapply(out, FUN = function(m) {
  ## whoops: 'lamb.int' vs 'lambda.int'
      colnames(m) <- gsub(colnames(m), pattern = "lamb\\.int\\[(.*)", replace = "lambda.int[\\1")
      nm.v <- sapply(strsplit(colnames(m), split = "\\["), FUN = "[", 1)
      as.mcmc(m[,nm.v %in% want])
    }))
  out.stack <- do.call(rbind, sub.out)
  obj.maker.dimAdd(out.stack)
},
  mc.cores = 4)   ## very memory-limited!

#
####
#################
################################################   BLOCK 2 -- end  #############



################################################   BLOCK 3 -- begin  ###########
#################    Calculate stats
####

nm.lut <- cbind(link = c("phi.int", "rho.int"), natural = c("phi.mn", "rho.mn"))

### empirical true values:
emp.lst <- mclapply(1:length(fl), FUN = function (n) {
  load(paste0(out.dir, "/", fl[n]))
  tmp <- sim.lst$emp
  names(tmp) <- nm.lut[match(names(names(tmp)), nm.lut[,"natural"]),"link"]    ## yikes, fragile !!!! names(names())
  tmp$phi.int <- qlogis(tmp$phi.int)
  tmp$rho.int <- log(tmp$rho.int)  
  tmp
})

emp.m.l <- sapply(names(emp.lst[[1]]), FUN = function (nm) {
  tmp <- lapply(emp.lst, FUN = get, x = nm)
  do.call(abind, args = tmp)
})


###  now abind() all these elements together:
###  this might not be RIGHT !!
smush <- sapply(names(fit.lst[[1]]), FUN = function(nm) {
  sub.lst <- lapply(fit.lst, FUN = "[[", nm)
  Reduce(function(m1, m2) {
   abind(m1, m2, along = length(dim(m1))-1)
   },
   sub.lst)
  })

diff.lst <- sapply(want, FUN = function (nm) {
  smush[[nm]] - true[[nm]]
  })

smush.trans <- c(lapply(smush[intersect(logit.par, names(smush))], FUN = "plogis"),
                 lapply(smush[intersect(log.par, names(smush))], FUN = "exp"))

diff.trans.lst <- sapply(want, FUN = function (nm) {
  smush.trans[[nm]] - true.trans[[nm]]
  })

diff.trans.prop <- sapply(want, FUN = function (nm) {
  (smush.trans[[nm]] - true.trans[[nm]])/true.trans[[nm]]
  })
  
diff.emp <- sapply(c("phi.int","rho.int"), FUN = function (nm) {
  sweep(smush[[nm]], MAR = 1:(length(dim(smush[[nm]]))), STATS = emp.m.l[[nm]], FUN = `-`)   ### to do recycling in an array!
})

bias.med <- sapply(diff.lst, FUN = function (a) {
  ## first, mean bias for each MCMC sample, then median over those samples:
#   stopifnot(length(dim(a)) %in% 0:3)
  if (length(dim(a)) %in% 0:1) { return(NA) }
  if (length(dim(a)) == 2) { a <- array(a, dim = c(1,dim(a))) }
  tmp <- apply(a, MAR = c(1,3), FUN = mean)
  apply(tmp, MAR = 1, FUN = median)
  })

bias.mn <- sapply(diff.lst, FUN = function (a) {
  ## equivalent to: first, mean bias for each MCMC sample, then mean over those samples:
  if (length(dim(a)) < 3) { a <- array(a, dim = c(1,dim(a))) }
  apply(a, MAR = 1, FUN = mean)
  })

rmse.med <- sapply(diff.lst, FUN = function (a) {
#   stopifnot(length(dim(a)) %in% 0:3)
  if (length(dim(a)) %in% 0:1) { return(NA) }
  if (length(dim(a)) == 2) { a <- array(a, dim = c(1,dim(a))) }
## first, mse for each MCMC sample, then median over those samples; then sqrt():
  tmp <- apply(a, MAR = c(1,3), FUN = function(v) { mean(v^2) })
  sqrt(apply(tmp, MAR = 1, FUN = median))
  })

rmse.mn <- sapply(diff.lst, FUN = function (a) {
#   stopifnot(length(dim(a)) %in% 0:3)
  if (length(dim(a)) %in% 0:1) { return(NA) }
  if (length(dim(a)) == 2) { a <- array(a, dim = c(1,dim(a))) }
## equivalent to: first, mse for each MCMC sample, then mean over those samples; then sqrt():
  sqrt(apply(a^2, MAR = c(1,3), FUN = mean))
  })

cvrg <- sapply(want, FUN = function (nm) {
  if (!is.null(smush[[nm]])) {
    lwr <- apply(smush[[nm]], MAR = 1:(length(dim(smush[[nm]]))-1), FUN = quantile, probs = c(0.025))
    upr <- apply(smush[[nm]], MAR = 1:(length(dim(smush[[nm]]))-1), FUN = quantile, probs = c(0.975))
  } else {
    lwr <- upr <- NA
  }
  
  tst <- true[[nm]] > lwr & true[[nm]] < upr
  if (is.null(dim(tst))) {
    sum(tst)
  } else {
    rowSums(tst)
  }
})



out.lst <- mclapply(1:length(fl), FUN = function (n) {
  load(paste0(out.dir, "/", fl[n]))
  as.mcmc.list(lapply(out, FUN = function(m) {
  ## whoops: 'lamb.int' vs 'lambda.int'
      colnames(m) <- gsub(colnames(m), pattern = "lamb\\.int\\[(.*)", replace = "lambda.int[\\1")
      nm.v <- sapply(strsplit(colnames(m), split = "\\["), FUN = "[", 1)
      as.mcmc(m[,nm.v %in% want])
    }))
  })
  
appnd.lst <- lapply(1:6, FUN = function (n) {
    as.mcmc(do.call(rbind, lapply(out.lst, FUN = `[[`, n)))
  })
  

mpsrf.lst <- gelman.diag(as.mcmc.list(appnd.lst))$mpsrf


#
####
#################
################################################   BLOCK 3 -- end  #############

return(mget(c("smush", "smush.trans", "diff.lst", "diff.trans.lst", "diff.trans.prop", "diff.emp", "bias.mn", "bias.med", "rmse.mn", "rmse.med", "cvrg", "true", "mpsrf.lst")))

}  ## END function body
