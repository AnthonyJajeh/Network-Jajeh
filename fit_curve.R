library(rjags)
library(coda)
library(tidyverse)

df <- read_csv("/Users/anthonyjajeh/Desktop/Research/Network-git/permeability_eps_data.csv")


# ensure plain numeric vectors
y   = as.numeric(log(df$k))   # log permeability
EPS = as.numeric(df$EPS)
n   =length(y)
lognormal_regression <- "
model {
  for(i in 1:n){
    y[i] ~ dnorm(mu[i], tau)         ## likelihood on log(k)
    mu[i] <- b[1] + b[2] * EPS[i]    ## linear predictor
  }

  ## priors
  for(j in 1:2){
    b[j] ~ dnorm(0, 1.0E-4)          ## weak priors for intercept & slope
  }
  tau ~ dgamma(0.001, 0.001)         ## residual precision
  sigma <- 1/sqrt(tau)

  ## derived parameters
  k0 <- exp(b[1])   ## permeability at EPS=0
  a  <- -b[2]       ## exponential decay rate
}
"

# load data

y = log(df$k)         # log permeability
EPS=df$EPS
n= length(y)
data = list(y=y, EPS=EPS, n=n)
inits <- list(b=c(0,0))

# run JAGS
j.model <- jags.model(file=textConnection(lognormal_regression),
                      data=data,
                      inits=inits,
                      n.chains=3)

update(j.model, 20000)   # burn-in

jags.out <- coda.samples(model=j.model,
                         variable.names=c("b","k0","a","sigma"),
                         n.iter=10000)

plot(jags.out)
print(gelman.diag(jags.out, multivariate = FALSE))
summary(jags.out)
samp <- sample.int(nrow(jags.out[[1]]), nsamp)

xpred = seq(min(EPS), max(EPS), length.out=50)
npred = length(xpred)
ycred = matrix(0, nsamp, npred)
ypred = matrix(0, nsamp, npred)

for(i in seq_len(nsamp)){
  b0 = jags.burn[[1]][samp[i], "b[1]"]
  b1 = jags.burn[[1]][samp[i], "b[2]"]
  sig= jags.burn[[1]][samp[i], "sigma"]
  
  mu = b0 + b1*xpred
  ycred[i,] = exp(mu)                         # mean curve (credible)
  ypred[i,] = exp(rnorm(npred, mu, sig))      # noisy curve (prediction)
}

ci = apply(ycred, 2, quantile, c(0.025,0.5,0.975))
pi = apply(ypred, 2, quantile, c(0.025,0.975))

plot(EPS, df$k, log="y", pch=16, col="darkgray")
lines(xpred, ci[2,], lwd=2, col="blue")  # median
lines(xpred, ci[1,], lty=2, col="blue")
lines(xpred, ci[3,], lty=2, col="blue")
lines(xpred, pi[1,], lty=2, col="red")
lines(xpred, pi[2,], lty=2, col="red")