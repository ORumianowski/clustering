library(nimble)
library(MASS)   # pour simuler des gaussiennes

# --------------------------------------------------
# 1. Données simulées
# --------------------------------------------------
set.seed(123)
N <- 200
D <- 2
K <- 2

true_mu <- list(c(-2, 0), c(3, 3))
true_sigma <- list(diag(2)*0.5, diag(2))
true_pi <- c(0.4, 0.6)

z <- rbinom(N, 1, true_pi[2]) + 1
X <- t(sapply(1:N, function(i) {
  if(z[i]==1) MASS::mvrnorm(1, true_mu[[1]], true_sigma[[1]])
  else MASS::mvrnorm(1, true_mu[[2]], true_sigma[[2]])
}))

# Fonction effort d'échantillonnage Ω(x)
effort_fun <- function(x) {
  # sous-échantillonnage si x[1] < 0
  if(x[1] < 0) return(0.3)
  return(1.0)
}
Omega <- apply(X, 1, effort_fun)

# --------------------------------------------------
# 2. Distribution personnalisée NIMBLE
# --------------------------------------------------
dGMM_weighted <- nimbleFunction(
  run = function(x = double(1),
                 mu = double(2),        # K x D
                 Sigma = double(3),     # D x D x K
                 pi = double(1),        # K
                 Omega = double(0),     # scalaire effort pour ce point
                 K = integer(0),
                 D = integer(0),
                 log = integer(0, default = 0)) {
    
    returnType(double(0))
    
    # densité brute du mélange
    dens <- 0.0
    for(k in 1:K) {
      dens <- dens + pi[k] * dmvnorm(x, mu[k,1:D], Sigma[1:D,1:D,k], log = FALSE)
    }
    
    # pondération locale
    dens_weighted <- Omega * dens
    
    # Approx Monte Carlo pour Z(theta,Omega)
    # Ici: on approxime sur une petite grille (peut être raffiné)
    M <- 50
    grid <- matrix(rnorm(M*D), M, D)  # points simulés ~ N(0, I) (simplification)
    dens_grid <- 0.0
    for(m in 1:M) {
      dtmp <- 0.0
      for(k in 1:K) {
        dtmp <- dtmp + pi[k] * dmvnorm(grid[m,], mu[k,1:D], Sigma[1:D,1:D,k], log = FALSE)
      }
      dens_grid <- dens_grid + effort_fun(grid[m,]) * dtmp
    }
    Z <- dens_grid / M
    
    val <- dens_weighted / (Z + 1e-12)
    if(log) return(log(val + 1e-12))
    return(val)
  }
)

# Enregistrement pour usage dans nimbleCode
registerDistributions(list(
  dGMM_weighted = list(
    BUGSdist = "dGMM_weighted(mu, Sigma, pi, Omega, K, D)",
    types = c("value = double(1)")
  )
))

# --------------------------------------------------
# 3. Modèle NIMBLE
# --------------------------------------------------
code <- nimbleCode({
  for(i in 1:N) {
    X[i,1:D] ~ dGMM_weighted(mu[1:K,1:D], Sigma[1:D,1:D,1:K], pi[1:K], Omega[i], K, D)
  }
  
  pi[1:K] ~ ddirch(alpha[1:K])
  
  for(k in 1:K) {
    mu[k,1:D] ~ dmnorm(mu0[1:D], cov = Sigma0[1:D,1:D])
    Sigma[1:D,1:D,k] ~ dwish(R[1:D,1:D], df)
  }
})

# --------------------------------------------------
# 4. Données et constantes
# --------------------------------------------------
constants <- list(N=N, D=D, K=K,
                  alpha=c(1,1),
                  mu0=c(0,0),
                  Sigma0=diag(D),
                  R=diag(D),
                  df=D+1)

data <- list(X = X, Omega = Omega)

inits <- list(
  mu = matrix(rnorm(K*D), K, D),
  Sigma = array(diag(D), dim = c(D,D,K)),
  pi = rep(1/K, K)
)

# --------------------------------------------------
# 5. Compilation et MCMC
# --------------------------------------------------
model <- nimbleModel(code, data=data, constants=constants, inits=inits)
cmodel <- compileNimble(model)

conf <- configureMCMC(model, monitors = c("mu","pi"))
mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc, project=model)

samples <- runMCMC(cmcmc, niter=2000, nburnin=500, thin=5)
print(head(samples))
