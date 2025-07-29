
data {
## Standard log-linear and logit-linear priors for normal distribution precision:  
  log.tau <- 0.15
  logit.tau <- 0.37
}

model {

## priors:
lambda.int ~ dnorm(0, log.tau)
a.int ~ dnorm(0, logit.tau)
r.int ~ dnorm(0, log.tau)
iota.int ~ dnorm(0, log.tau)
## p.int is coming in as data!

log(lambda) <- lambda.int
logit(a) <- a.int
log(r) <- r.int
log(iota) <- iota.int
logit(p) <- p.int

for (s in 1:n.site) {
  N[s,1] ~ dpois(lambda)
  for (t in 2:n.yr) {
    S[s,t] ~ dbinom(a, N[s,t-1])
    G[s,t] ~ dpois(N[s,t-1] * r)
    I[s,t] ~ dpois(iota)
    N[s,t] <- S[s,t] + G[s,t] + I[s,t]
    y[s,t] ~ dbinom(p, N[s,t])  
  }  ## END t loop
}  ## END s loop

}