##  Collect and plot results of fitting model to a number of simulated datasets.
##
##
##
############################################

## loads function gathrFits():
source('functions/gather_simFits.r')

library(vioplot)


## Output directory:
#out.dir <- 'output/mcmc'  ## the default

key.m <- cbind(type = "effBreed",
               nm = c("", "brPhi_", "brPhiNbPhi", "brPhiRho", "brNbPhiRho_2", "brNbPhiRho_tight"),
               mod.date = c("05dec", "06dec", "23dec", "06dec", "23dec", "23dec"),
               fit.date1 = c("12-22", "12-22", "12-2[3|4]", "12-23", "12-26", "12-26"))


rslt <- apply(key.m, MAR = 1, FUN = gathrFits)
names(rslt) <- key.m[,"nm"]

## mod letters:
key.m <- cbind(key.m, mod.lett = LETTERS[1:nrow(key.m)])

## convergence:
mpsrf.v <- sapply(rslt, FUN = function (l) { l$mpsrf.lst })
names(mpsrf.v) <- LETTERS[1:nrow(key.m)]


### Make a df with 'diff' values, and label of model: no 'Psi'
toPlot <- data.frame(par = c(rep("lambda.int",2), rep("phi.int", 4), rep("eps.int", 2), rep("iota.int", 2), "rho.int"),
                addr = c(1:2, 1:4, 1:2, 1:2, 1)) 

diff.lst <- lapply(rslt, FUN = `[[`, "diff.lst")
diff.trans.lst <- lapply(rslt, FUN = `[[`, "diff.trans.lst")
diff.trans.prop <- lapply(rslt, FUN = `[[`, "diff.trans.prop")
smush.trans <- lapply(rslt, FUN = `[[`, "smush.trans")

names(diff.lst) <- names(diff.trans.lst) <- names(diff.trans.prop) <- names(smush.trans) <- key.m[,"mod.lett"]

### Should all be the same... right?
dm <- dim(diff.lst[[1]]$rho.int)

m <- apply(toPlot, MAR = 1, FUN = function (v) {   ### transformed, through link function
  sapply(diff.lst, FUN = function (l) {
    a <- l[[as.character(v["par"])]]
    if( length(dim(a)) == 2 ) {
      a <- array(a, dim = c(1,dim(a)))
    } 
    a[as.integer(v["addr"]),,]
  })})

m.trans <- apply(toPlot, MAR = 1, FUN = function (v) {  ## confusingly, this is in fact back-transformed
  sapply(diff.trans.lst, FUN = function (l) {
    a <- l[[as.character(v["par"])]]
    if( length(dim(a)) == 2 ) {
      a <- array(a, dim = c(1,dim(a)))
    } 
    a[as.integer(v["addr"]),,]
  })})

m.prop <- apply(toPlot, MAR = 1, FUN = function (v) {   ## back-transformed, then difference expressed as a proportion of true value
  sapply(diff.trans.prop, FUN = function (l) {
    a <- l[[as.character(v["par"])]]
    if( length(dim(a)) == 2 ) {
      a <- array(a, dim = c(1,dim(a)))
    } 
    a[as.integer(v["addr"]),,]
  })})


four.seas <- c("br", "m1", "nb", "m2")
two.seas <- four.seas[c(1,3)]

par.nm <- cbind(toPlot,
                full = c(rep("lambda", 2), rep("phi", 4), rep("epsilon", 2), rep("iota", 2), "rho"),
                seas = c(paste0("[", c(two.seas, four.seas, two.seas, two.seas), "]"),  ""))

colnames(m.prop) <- colnames(m.trans) <- colnames(m) <- apply(par.nm[c("full","seas")], MAR = 1, FUN = function (v) {
  paste0(v[1], v[2])
  })


## or:
rownames(m.prop) <- rownames(m.trans) <- rownames(m) <- rep(names(diff.lst), each = prod(dm))


## alternatively, could put all parameters together in the same panel, for each model:
vio.fun <- function(mat, mod, par, col.v, xlimit, max.vals = 1e5, ylabs = TRUE, yside = 2, y.cex = 1.5, med.ramp = FALSE) {
  cols <- rbind(c(0.5,0.5,0.8),
                c(0.8,0.5,0.5))
  sub <- mat[rownames(mat) == mod, par, drop=FALSE][1:max.vals,]
  df <- data.frame(par = rep(colnames(sub), each = nrow(sub)), val = c(sub))
  vioplot(val ~ factor(par, levels = rev(colnames(sub))), data = df, col = "white", border = "white", rectCol = "white", lineCol = "white", horizontal = TRUE, xaxt = "n", ylab = "", xlab = "", ylim = xlimit)
  if (ylabs) { axis(side = yside, labels = parse(text=rev(colnames(sub))), at = 1:length(colnames(sub)), las = 1, cex.axis = y.cex) }
  abline(v=0)
  vp <- vioplot(val ~ factor(par, levels = rev(colnames(sub))), data = df,
    col = rgb(cols[col.v,], alpha = 0.6),
    border = rgb(cols[col.v,], alpha = 0.8),
    rectCol = rgb(0.5,0.5,0.5),  ## alpha doesn't seem to work in this feature
    lineCol = rgb(0.5,0.5,0.5),  ## alpha doesn't seem to work in this feature
    horizontal = TRUE, add = T)
  if (med.ramp) { points(x = vp$median, y = 1:length(vp$median), pch = 16, col = rgb({tmp <- rmp(abs(vp$median)/xlimit[2]); tmp[is.na(tmp)] <- 0; tmp}, maxC = 255)) }
}

rn <- colnames(m)          #s# parameters
cn <- unique(rownames(m))  ## models

prior.m <- array(1, dim = lengths(list(rn,cn)), dimnames = list(par = rn, mod = cn))

## make these '2':
prior.m["phi[br]",2:6] <- 2              # br phi
prior.m["rho",4:6] <- 2             # rho
prior.m["phi[nb]",c(3,5:6)] <- 2              # nb phi
#prior.m["epsilon[br]",7] <- 2

##  color ramp:
### get good ColorBrewer colors to make a ramp from:
library(RColorBrewer)
ends <- brewer.pal(5, "Greens")[c(1,5)]
rmp <- colorRamp(ends, space = "Lab")

nms <- unique(rownames(m))

par.log <- c("lambda[br]", "lambda[nb]", "iota[br]", "iota[nb]", "rho")
par.logit <- c("phi[br]", "phi[m1]", "phi[nb]", "phi[m2]", "epsilon[br]", "epsilon[nb]")


### one on left side, one on right:
lab.df <- data.frame(lab = c(TRUE, rep(FALSE, length(nms)-2), TRUE),
                      side = c(2, rep(NA, length(nms)-2), 4))

### all on left side:
lab.df <- data.frame(lab = rep(TRUE, length(nms)),
                      side = rep(2, length(nms)))

##########################################################################################
#############################################
##  FIGURE:

png(file = paste0('figs/violin/violin_', Sys.Date(), '_twoRow.png'), width = 15, height = 7, units = "in", res = 400)
### for axis labels on right side:
#par(mfcol = c(2, length(nms)), mar = c(3,3,1,1), oma = c(2,3,0,5))
### for NO axis labels on right side:
par(mfcol = c(2, length(nms)), mar = c(3,3,1,1), oma = c(2,3,0,0))

sapply(1:length(nms), FUN = function (s) {

## logit:
  vio.fun(mat = m, mod = nms[s], par = par.logit, col.v = rev(prior.m[par.logit,nms[s]]), xlimit = c(-3,3), ylabs = lab.df[s,"lab"], yside = lab.df[s,"side"], y.cex = 1.5)
  u <- par("usr")
  text(x = u[1] + 0.9*(u[2]-u[1]), y = u[3] + 0.95*(u[4]-u[3]), label = nms[s], cex = 1.8)

## log:  
  vio.fun(mat = m, mod = nms[s], par = par.log[grep(par.log, pattern = "lambda", invert = TRUE)], col.v = rev(prior.m[par.log,nms[s]]), xlimit = c(-1,1), ylabs = lab.df[s,"lab"], yside = lab.df[s,"side"], y.cex = 1.5)
  u <- par("usr")
  text(x = u[1] + 0.9*(u[2]-u[1]), y = u[3] + 0.95*(u[4]-u[3]), label = nms[s], cex = 1.8)

})

mtext(side = 1, outer = TRUE, line = 1, text = "Difference from true value")
mtext(side = 2, outer = TRUE, line = 1.5, text = "Transformed parameter")
dev.off()


##########################################################################################
#############################################
##  TABLE (Supplement):


##  Get median bias, median rmse, coverage:
bias.med <- lapply(rslt, FUN = `[[`, "bias.med")
rmse.med <- lapply(rslt, FUN = `[[`, "rmse.med")
cvrg <- lapply(rslt, FUN = `[[`, "cvrg")

names(bias.med) <- names(rmse.med) <- names(cvrg) <- key.m[,"mod.lett"]

##  label list elements; will make things easier later:
par <- c("lambda.int", "phi.int", "eps.int", "iota.int")
seas <- list(c("br", "nb"), c("br", "m1", "nb", "m2"), c("br", "nb"), c("br", "nb"))

lab.m <- cbind(par = rep(par, lengths(seas)), seas = unlist(seas)


##  Arrange dataframe:
##  Make the parameters the rows, group estimates by scenario across the top:
df.l <- lapply(names(bias.med), FUN = function (l) {
   b.l <- bias.med[[l]]
   r.l <- rmse.med[[l]]
   c.l <- cvrg[[l]]

   tmp <- cbind(round(unlist(b.l[names(b.l) %in% want]), 2),
                round(unlist(r.l[names(r.l) %in% want]), 2),
                unlist(c.l[names(c.l) %in% want]))
})
names(df.l) <- names(bias.med)

meas <- do.call(cbind, df.l)
colnames(meas) <- paste(rep(names(df.l), each = 3), rep(c("Bias", "RMSE", "Coverage"), length(df.l)), sep = ".")

df <- data.frame(lab.m, meas)


library(xtable)

tab <- xtable(df)

print(tab, include.rownames = FALSE)


% latex table generated in R 4.3.3 by xtable 1.8-4 package
% Thu Jun  6 18:04:53 2024
\begin{table}[ht]
\centering
\begin{tabular}{llrrrrrrrrrrrrrrrrrr}
  \hline
par & seas & A.Bias & A.RMSE & A.Coverage & B.Bias & B.RMSE & B.Coverage & C.Bias & C.RMSE & C.Coverage & D.Bias & D.RMSE & D.Coverage & E.Bias & E.RMSE & E.Coverage & F.Bias & F.RMSE & F.Coverage \\
  \hline
lambda.int & br & -0.00 & 0.07 & 98 & -0.01 & 0.08 & 95 & -0.01 & 0.08 & 96 & 0.00 & 0.07 & 98 & 0.00 & 0.07 & 97 & -0.01 & 0.07 & 98 \\
  lambda.int & nb & -0.01 & 0.17 & 95 & -0.02 & 0.18 & 93 & 0.00 & 0.17 & 96 & -0.02 & 0.17 & 94 & -0.00 & 0.17 & 94 & -0.01 & 0.17 & 93 \\
  phi.int & br & -1.98 & 2.53 & 100 & 0.52 & 1.41 & 100 & 0.53 & 1.42 & 100 & 0.44 & 1.40 & 100 & 0.46 & 1.39 & 100 & 0.15 & 0.67 & 100 \\
  phi.int & m1 & -0.49 & 1.51 & 100 & -0.63 & 1.60 & 100 & -0.70 & 1.69 & 100 & 0.40 & 1.10 & 100 & 0.37 & 1.08 & 100 & 0.39 & 1.05 & 100 \\
  phi.int & nb & -0.55 & 1.56 & 100 & -0.69 & 1.63 & 100 & 0.19 & 0.81 & 100 & 0.41 & 1.11 & 100 & 0.36 & 0.82 & 100 & 0.10 & 0.42 & 100 \\
  phi.int & m2 & -0.88 & 1.71 & 100 & -1.01 & 1.81 & 100 & -0.97 & 1.80 & 100 & 0.13 & 1.04 & 100 & 0.07 & 1.00 & 100 & 0.10 & 0.98 & 100 \\
  eps.int & br & 0.14 & 1.26 & 100 & 0.21 & 1.28 & 100 & 0.29 & 1.34 & 100 & -0.34 & 1.09 & 100 & -0.32 & 1.07 & 100 & -0.34 & 1.05 & 100 \\
  eps.int & nb & -0.17 & 0.95 & 100 & -0.05 & 0.94 & 100 & -0.15 & 0.95 & 100 & -0.05 & 0.89 & 97 & -0.12 & 0.88 & 98 & 0.02 & 0.86 & 100 \\
  iota.int & br & -2.02 & 2.37 & 38 & -1.53 & 1.90 & 50 & -1.12 & 1.56 & 69 & 0.01 & 0.16 & 99 & -0.01 & 0.16 & 100 & -0.00 & 0.16 & 97 \\
  iota.int & nb & -1.81 & 2.14 & 34 & -1.35 & 1.70 & 43 & -1.03 & 1.42 & 61 & -0.05 & 0.11 & 98 & -0.03 & 0.11 & 98 & -0.02 & 0.10 & 100 \\
   \hline
\end{tabular}
\end{table}
