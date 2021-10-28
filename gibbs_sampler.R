# Libraries ---------------------------------------------------------------

library(mvnfast)

# Data Simulation ---------------------------------------------------------

source("data_simulation.R")

# Initialization ----------------------------------------------------------

al <- 0.5 * rep(1, km) # gamma hyperparameters for latent variable precisions
bl <- 0.5 * rep(1, km) # gamma hyperparameters for latent variable precisions
cj <- 1 * rep(1, p)    # gamma hyperparameters for error precisions
dj <- 0.2 * rep(1, p)  # gamma hyperparameters for error precisions

pl <- km * (p - km) + km * (km + 1) / 2  # number of free parameters in factor matrix
nc <- c(1:km, rep(km, p - km)) # number of free parameters in each row of Lambda

Plam <- diag(pl)         # used to specify the prior precision matrix for Lambda^*
Pe   <- diag(p)          # initial value for error precision matrix
Pl   <- 3 * diag(km)     # initial precision matrix for latent variables
Ls   <- matrix(0, p, km) # initial value for Lambda^*
Ls[1:km, 1:km] <- diag(km)

# MCMC parameters and arrays ----------------------------------------------

n_sim     <- 25000 # total number of MCMC iterations
n_burn_in <- 5000  # burn-in for Gibbs

# Define arrays for saving output
Lout <- matrix(0, n_sim - n_burn_in, p * km) # for Lambda
Pout <- matrix(0, n_sim - n_burn_in, p)      # for Sigma^(-1)
Oout <- matrix(0, n_sim - n_burn_in, p * p)  # for Omega

# Sample from full conditionals -------------------------------------------

for (iter in 1:n_sim) {
  
  # Step 1 - update latent factors
  
  # covariance matrix for factors
  Veta <- solve(Pl + t(Ls) %*% Pe %*% Ls)   
  
  # mean vector for factors
  Eeta     <- t(Veta %*% t(Ls) %*% Pe %*% t(Y)) 
  eta      <- apply(Eeta, 1, function(x) rmvn(1, x, Veta))
  dim(eta) <- c(n, km)
  # eta  <- rmvn(nrow(Eeta), Eeta, Veta)      # latent factors
  
  # Step 2 - update factor loadings
  
  for (j in 1:p) { # jth row of Lambda
    
    # etaS <- eta[, 1:nc[j]]
    z_j      <- eta[, 1:nc[j]]
    dim(z_j) <- c(n, min(j, km))
    
    # covariance matrix for loadings
    Vlam <- solve(
      Plam[nc[j], nc[j]] + Pe[j, j] * t(z_j) %*% z_j
    ) 
    
    # mean vector for for loadings
    Elam <- Vlam %*% (Pe[j, j] * t(z_j) %*% Y[, j, drop = FALSE])
    
    # Lambda^*: factor loadings under PX model
    Ls[j, 1:nc[j]] <- rmvn(1, Elam, Vlam)
  }
  
  # Step 3 - update latent variable precision
  ae <- al + n / 2
  be <- bl + 0.5 * t(eta ^ 2) %*% matrix(1, n, 1)
  
  # latent variable precision matrix
  Pl <- diag(rgamma(km, ae, be), nrow = km, ncol = km) 
  
  # Step 4 - update residual precision
  ap <- cj + n / 2
  bp <- dj + 0.5 * t((Y - eta %*% t(Ls)) ^ 2) %*% matrix(1, n, 1)
  
  # error precision matrix
  Pe <- diag(rgamma(p, ap, bp), nrow = p, ncol = p) 
  
  # Step 5 - Recalculate original factor loadings and save sampled values
  L <- Ls
  
  for (j in 1:km) {
    if (Ls[j, j] < 0) {
      L[, j] <- -L[, j]
    }
  }
  
  L <- L %*% sqrt(solve(Pl))
  
  if (iter > n_burn_in) {
    Lout[iter - n_burn_in, ] <- c(L)
    Pout[iter - n_burn_in, ] <- diag(Pe)
    Oout[iter - n_burn_in, ] <- c(L %*% t(L) + solve(Pe))
  }
  
  if (iter %% 1000 == 0 & iter <= n_burn_in) {
    cat(paste0("Iteration: ", iter, " (burn in).\n"))
  } else if (iter %% 1000 == 0 & iter > n_burn_in) {
    cat(paste0("Iteration: ", iter, " (sampling).\n"))
  }
}

readr::write_rds(
  list(
    Lout = Lout,
    Pout = Pout,
    Oout = Oout
  ),
  "data/posterior_samples.rds"
)
