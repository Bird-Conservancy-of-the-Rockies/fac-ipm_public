##  Collect convergence diagnostics of fitting model to a number of simulated datasets.
##
##
##  14 dec 2023
############################################

gathrConv <- function (
  key,       ## vector of keywords to look for in files
  out.dir = 'output/mcmc',
  max.n.files = 10
) {

require(coda)
require(parallel)
require(abind)

#source('../fit/functions/obj_maker_dimAdd.r')



lf <- list.files(out.dir)
tst.m <- sapply(key, FUN = function (s) { grepl(lf, pattern = s, fixed = T) })
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

#logit.par <- c("psi.int", "phi.int", "eps.int")
#log.par <- c("lambda.int", "iota.int", "rho.int")

#true.trans <- c(lapply(true[logit.par], FUN = plogis),
#                lapply(true[log.par], FUN = exp))


#
####
#################
################################################   BLOCK 1 -- end  #############


################################################   BLOCK 2 -- begin  ###########
#################    Collect all MCMC results
####
#  In an array?  I guess I want to know the bias / coverage for each dataset...

want <- names(true)

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

  
#  sub.out <- as.mcmc.list(lapply(out, FUN = function(m) {
#  ## whoops: 'lamb.int' vs 'lambda.int'
#      colnames(m) <- gsub(colnames(m), pattern = "lamb\\.int\\[(.*)", replace = "lambda.int[\\1")
#      nm.v <- sapply(strsplit(colnames(m), split = "\\["), FUN = "[", 1)
#      as.mcmc(m[,nm.v %in% want])
#    }))
##  g <- gelman.diag(sub.out)
#  out.stack <- do.call(rbind, sub.out)
##  obj.maker.dimAdd(out.stack)
#  as.mcmc(out.stack)
#},
#  mc.cores = 4)   ## very memory-limited!

#gelman.diag(as.mcmc.list(mcmc.lst))
#mcmc.stack <- do.call(rbind, mcmc.lst)


#conv.psrf <- do.call(abind, args = c(lapply(conv.lst, FUN = function(l) { l$psrf }), along = 3) )

#conv.mpsrf <- sapply(conv.lst, FUN = function(l) { l$mpsrf }) 


return(mget("mpsrf.lst"))

}  ## END function body
