

# produit des erreurs mais sinon le corps y est

library(nimble)
library(MASS)

# ----- PARAMÈTRES -----
set.seed(123)
N <- 150
K <- 3
d <- 4
Tmax <- 5

# Proportions et clusters
pi_true <- c(0.3, 0.4, 0.3)
z_true <- sample(1:K, N, replace=TRUE, prob=pi_true)

# Moyennes et variances
mu_true <- list(c(0,0,0,0), c(3,3,3,3), c(-3,3,-3,3))
Sigma_true <- list(diag(1,d), diag(1,d), diag(1,d))

# Survie (constante dans le temps) et détection
phi_true <- c(0.9, 0.7, 0.5)
p_true <- 0.6

# Coordonnées
X <- t(sapply(1:N, function(i) mvrnorm(1, mu_true[[z_true[i]]], Sigma_true[[z_true[i]]])))

# États latents et observations
state <- matrix(NA, nrow=N, ncol=T)  # 1 = vivant, 2 = mort
y <- matrix(NA, nrow=N, ncol=T)      # 1 = détecté, 2 = non détecté

for(i in 1:N){
  k <- z_true[i]
  state[i,1] <- 1  # tous vivants au départ
  for(t in 2:Tmax){
    state[i,t] <- sample(1:2, size=1, prob=c(phi_true[k], 1 - phi_true[k])) * (state[i,t-1] == 1) + 2 * (state[i,t-1] == 2)
  }
  for(t in 1:Tmax){
    if(state[i,t] == 1){
      y[i,t] <- sample(1:2, size=1, prob=c(p_true, 1 - p_true))
    } else {
      y[i,t] <- 2  # jamais détecté si mort
    }
  }
}

# ----- CODE NIMBLE -----

code <- nimbleCode({
  pi[1:K] ~ ddirch(alpha[1:K])
  p ~ dbeta(1,1)
  
  for(k in 1:K){
    mu[k,1:d] ~ dmnorm(mean=mu0[1:d], prec=Prec0[1:d,1:d])
    Prec[k,1:d,1:d] ~ dwish(R[1:d,1:d], df)
    Sigma[k,1:d,1:d] <- inverse(Prec[k,1:d,1:d])
    
    phi[k] ~ dbeta(1,1)  # survie constante par cluster
    
    # Matrice de transition 2x2 pour chaque cluster
    transition[k,1,1] <- phi[k]
    transition[k,1,2] <- 1 - phi[k]
    transition[k,2,1] <- 0
    transition[k,2,2] <- 1
  }
  
  # Matrice d'observation (commune à tous)
  obs[1,1] <- p       # vivant -> détecté
  obs[1,2] <- 1 - p   # vivant -> non détecté
  obs[2,1] <- 0       # mort -> détecté
  obs[2,2] <- 1       # mort -> non détecté
  
  for(i in 1:N){
    z[i] ~ dcat(pi[1:K])
    
    # Recalcul des moyennes et précisions (via variable intermédiaire)
    for(j in 1:d){
      mu_i[i,j] <- mu[z[i], j]
      for(k1 in 1:d){
        Prec_i[i,j,k1] <- Prec[z[i], j, k1]
      }
    }
    
    X[i,1:d] ~ dmnorm(mu_i[i,1:d], prec=Prec_i[i,1:d,1:d])
    
    state[i,1] ~ dcat(init[1:2])  # tous vivants au départ
    
    for(t in 2:Tmax){
      state[i,t] ~ dcat(transition[z[i], state[i,t-1], 1:2])
    }
    
    for(t in 1:Tmax){
      y[i,t] ~ dcat(obs[state[i,t], 1:2])
    }
  }
})

# ----- INITIALISATION -----

constants <- list(
  N = N, K = K, d = d, Tmax = Tmax,
  alpha = rep(1,K),
  R = diag(d), df = d+1,
  mu0 = rep(0,d), Prec0 = diag(0.01,d),
  init = c(1, 0)  # vivant au temps 1
)

data <- list(X = X, y = y)

# inits <- list(
#   z = sample(1:K, N, replace=TRUE),
#   mu = matrix(0, K, d),
#   Prec = array(rep(diag(d), K), dim=c(K,d,d)),
#   phi = rep(0.5, K),
#   p = 0.5,
#   state = matrix(1, N, Tmax)
# )


inits <- list(
  z = z_true,  # vrais clusters
  mu = do.call(rbind, mu_true),  # vraies moyennes
  Prec = array(unlist(lapply(Sigma_true, solve)), dim = c(K, d, d)),  # précision = inverse de Sigma
  phi = phi_true,
  p = p_true,
  state = state  # états latents vrais
)
# ----- MCMC -----

model <- nimbleModel(code, data=data, constants=constants, inits=inits)
Cmodel <- compileNimble(model)

conf <- configureMCMC(model, monitors=c("pi", "mu", "Sigma", "phi", "p"))
Rmcmc <- buildMCMC(conf)
Cmcmc <- compileNimble(Rmcmc, project=model)

samples <- runMCMC(Cmcmc, niter=3000, nburnin=1000, thin=5, summary=TRUE)

# Résultats
print(samples$summary)
