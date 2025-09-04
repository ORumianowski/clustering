
library(nimble)


# Simulation des données --------------------------------------------------

set.seed(123)
N <- 150
K <- 3
d <- 4

# Proportions et clusters
pi_true <- c(0.3, 0.4, 0.3)
z_true <- sample(1:K, N, replace=TRUE, prob=pi_true)

# Moyennes et covariances
mu_true <- list(c(0,0,0,0), c(3,3,3,3), c(-3,3,-3,3))
Sigma_true <- list(diag(1,d), diag(1,d), diag(1,d))

# Survie par cluster
phi_true <- c(0.9, 0.7, 0.5)

# Simuler les coordonnées
X <- t(sapply(1:N, function(i) MASS::mvrnorm(1, mu_true[[z_true[i]]], Sigma_true[[z_true[i]]])))

# Simuler la survie
y <- sapply(1:N, function(i) rbinom(1,1,phi_true[z_true[i]]))



# Modèle hiérarchique latent (CMR avec clusters) --------------------------

# =======================
# Modèle intégré clusters 4D + survie (corrected)
# =======================
code <- nimbleCode({
  # Prior sur les proportions de clusters
  pi[1:K] ~ ddirch(alpha[1:K])
  
  # Paramètres gaussiens des clusters
  for(k in 1:K){
    mu[k,1:d] ~ dmnorm(mean = rep(0,d), prec = diag(0.01,d))  # prior vague
    Prec[k,1:d,1:d] ~ dwish(R[1:d,1:d], df)
    Sigma[k,1:d,1:d] <- inverse(Prec[k,1:d,1:d])
    
    # Survie par cluster
    phi[k] ~ dbeta(1,1)
  }
  
  # Latent cluster et observations
  for(i in 1:N){
    z[i] ~ dcat(pi[1:K])
    
    # Définir mu_i et Prec_i comme vecteurs et matrices déterministes
    for(j in 1:d){
      mu_i[i,j] <- mu[z[i],j]
    }
    for(j in 1:d){
      for(k1 in 1:d){
        Prec_i[i,j,k1] <- Prec[z[i],j,k1]
      }
    }
    
    # Observations
    X[i,1:d] ~ dmnorm(mean = mu_i[i,1:d], prec = Prec_i[i,1:d,1:d])
    y[i] ~ dbern(phi[z[i]])
  }
})




# Initialisation et données -----------------------------------------------

constants <- list(N=N, K=K, d=d, alpha=rep(1,K), R=diag(d), df=d+1)
data <- list(X=X, y=y)
inits <- list(
  z = sample(1:K, N, replace=TRUE),
  mu = matrix(0, nrow=K, ncol=d),
  Prec = array(rep(diag(d),K), dim=c(K,d,d)),
  phi = rep(0.5,K)
)



# Compilation et MCMC -----------------------------------------------------

model <- nimbleModel(code, data=data, constants=constants, inits=inits)
Cmodel <- compileNimble(model)

conf <- configureMCMC(model, monitors=c("pi","mu","Sigma","phi","z"))
Rmcmc <- buildMCMC(conf)
Cmcmc <- compileNimble(Rmcmc, project=model)

samples <- runMCMC(Cmcmc, niter=5000, nburnin=1000, thin=5)

