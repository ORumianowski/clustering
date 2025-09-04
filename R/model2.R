library(nimble)
library(MASS)
library(clue)

set.seed(123)

# ============================
# 1. Simulation des données
# ============================
K <- 3
n <- 300
d <- 4

true_means <- list(c(2,2,0,1),
                   c(3,0,2,1),
                   c(0,3,3,2))

true_covs <- list(diag(d),diag(d),diag(d))


true_probs <- c(0.3,0.4,0.3)
z_true <- sample(1:K, n, replace=TRUE, prob=true_probs)

X <- matrix(NA, nrow=n, ncol=d)
for(i in 1:n){
  X[i,] <- mvrnorm(1, mu=true_means[[z_true[i]]], Sigma=true_covs[[z_true[i]]])
}


# ============================
# 2. Définition du modèle NIMBLE
# ============================
code <- nimbleCode({
  pi[1:K] ~ ddirch(alpha[1:K])
  for(k in 1:K){
    mu[k,1:d] ~ dmnorm(mean=mu0[1:d], cov=cov0[1:d,1:d])
    Prec[k,1:d,1:d] ~ dwish(R[1:d,1:d], df)
  }
  for(i in 1:N){
    z[i] ~ dcat(pi[1:K])
    X[i,1:d] ~ dmnorm(mean=mu[z[i],1:d], prec=Prec[z[i],1:d,1:d])
  }
})

# ============================
# 3. Constantes et données
# ============================
constants <- list(K=K, N=n, d=d,
                  mu0=rep(0,d), cov0=diag(d)*10,
                  R=diag(d), df=d+1, alpha=rep(1,K))
data <- list(X=X)

km <- kmeans(X, centers=K)

inits <- list(
  z = km$cluster,
  mu = km$centers,
  Prec = array(rep(diag(d), K), dim=c(K,d,d))
)


# ============================
# 4. Compilation et MCMC
# ============================
model <- nimbleModel(code, data=data, constants=constants, inits=inits)
Cmodel <- compileNimble(model)

conf <- configureMCMC(model, monitors=c("pi","mu","Prec","z"))
Rmcmc <- buildMCMC(conf)
Cmcmc <- compileNimble(Rmcmc, project=model)

samples <- runMCMC(Cmcmc, niter=6000, nburnin=1000, thin=5, nchains=1, setSeed=123)
