##  A function that ingests an ENTIRE matrix of MCMC samples (with JAGS-given column names; i.e., an mcmc.list object will have been previously `stacked` as via `do.call(rbind, {mcmc.lst})`), and returns a list of objects each with its proper dimensions, plus one for either `n.iter` mcmc samples, or all the samples if `n.iter = NULL`. 
##
##  19 may 2022, BLN
##########################################


##  Intended to be source()'d: 
##  All it does is load a single function:

obj.maker.dimAdd <- function (mcmc, n.iter = NULL) {
  if (is.null(n.iter))  n.samp <- nrow(mcmc)
  mcmc.nms <- colnames(mcmc)
  nms <- gsub(mcmc.nms, pattern = "(*.)\\[.*", replace = "\\1")
## determine which rows came from a vector, which from a matrix (don't think there'll be 3d arrays in any nodes here)
  bracks <- gsub(mcmc.nms, pattern = ".*\\[(.*)\\]", replace = "\\1")
  two.comma <- grep(bracks, pattern = ",.+,")
  one.comma <- setdiff(grep(bracks, pattern = ","), two.comma)
  unq <- unique(nms)
  ind.lst <- sapply(unq, FUN = function (s) { which(nms == s) }, simplify = FALSE)
  three.d <- unique(nms[two.comma])
  two.d <- unique(nms[one.comma])
  scalar <- names(ind.lst)[sapply(ind.lst, FUN = function (v) { length(v) == 1 })]
  one.d <- unique(nms[!nms %in% c(three.d, two.d, scalar)])

##  Treat 3-d objects separately from 1-d and 2-d:
  xyz.lst <- sapply(three.d, FUN = function(s) { 
              ind <- ind.lst[[s]]
              split <- lapply(strsplit(bracks[ind], split = ","), FUN = as.integer)
              m <- do.call(rbind, split)
              cbind(ind = ind, x = m[,1], y = m[,2], z = m[,3])
            }, simplify = FALSE)

##  Treat 2-d objects separately from 1-d and 3-d:
  xy.lst <- sapply(two.d, FUN = function(s) { 
              ind <- ind.lst[[s]]
              split <- lapply(strsplit(bracks[ind], split = ","), FUN = as.integer)
              m <- do.call(rbind, split)
              cbind(ind = ind, x = m[,1], y = m[,2])
            }, simplify = FALSE)

##  One-dimensional:
  x.lst <- sapply(one.d, FUN = function (s) {
              ind <- ind.lst[[s]]              
              cbind(ind = ind, x = as.integer(bracks[ind]))
           }, simplify = FALSE)

##  Make the objects:
  xyz.obj <- sapply(xyz.lst, FUN = function (m) {
              arr <- array(NA, dim = c(max(m[,"x"]), max(m[,"y"]), max(m[,"z"]), n.samp), dimnames = list(x = NULL, y = NULL, z = NULL, iter = NULL))
              invisible(apply(m, MAR = 1, FUN = function (v) {
                      arr[v[2],v[3],v[4],] <<- mcmc[,v[1]] }))
              arr
            }, simplify = F)

  xy.obj <- sapply(xy.lst, FUN = function (m) {
              arr <- array(NA, dim = c(max(m[,"x"]), max(m[,"y"]), n.samp), dimnames = list(x = NULL, y = NULL, iter = NULL))
              invisible(apply(m, MAR = 1, FUN = function (v) {
                      arr[v[2],v[3],] <<- mcmc[,v[1]] }))
              arr
            }, simplify = F)

  x.obj <- sapply(x.lst, FUN = function (m) {
              mat <- matrix(NA, nrow = max(m[,"x"]), ncol = n.samp)
              invisible(apply(m, MAR = 1, FUN = function (v) {
                      mat[v[2],] <<- mcmc[,v[1]] }))
              mat
            }, simplify = F)
              
##  Scalars:
  sclr.obj <- sapply(scalar, FUN = function (s) {
                mcmc[,ind.lst[[s]]]
              }, simplify = F)


  return(c(sclr.obj, x.obj, xy.obj, xyz.obj, n.samp = n.samp))
}

