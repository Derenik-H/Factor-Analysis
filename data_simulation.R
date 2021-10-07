# Libraries ---------------------------------------------------------------

library(mvnfast)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

# Data simulation ---------------------------------------------------------

p     <- 7   # number of outcomes
n     <- 100 # sample size
truek <- 1   # true number of factors
km    <- 1   # maximum number of factors considered

Lt <- matrix(0, p, truek) # True Lambda

Lt[, 1] <- c(0.995, 0.975, 0.949, 0.922, 0.894, 0.866, 0.837)
St <- diag(c(0.01, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30)) # True Sigma

Ot <- Lt %*% t(Lt) + St # True Omega
mu <- rep(0, p) # True mu

Y <- rmvn(n, mu, Ot) # generate outcome variables, Y
Y <- scale(Y)        # centering and scaling the data
pairs(Y)
