# ============================================================
# Script : Clustering bayésien avec un modèle de mélange gaussien (NIMBLE)
# Objectif : Attribuer des points en 4 dimensions à K clusters
#            et estimer les paramètres (moyenne, covariance) des gaussiennes.
# ============================================================

# Installer les packages si besoin :
# install.packages("nimble")
# install.packages("MASS")

library(nimble)
library(MASS)

set.seed(123)

# ============================================================
# 1. Génération de données simulées (pour tester le modèle)
# ============================================================

K <- 3      # nombre de clusters fixé
n <- 300    # nombre de points
d <- 4      # dimension des données

# Paramètres "vrais" pour générer les données
true_means <- list(c(-2, -2, 0, 1),
                   c(3, 0, -2, -1),
                   c(0, 3, 3, -2))

true_covs <- list(diag(d),
                  matrix(c(2,0.5,0.3,0,
                           0.5,1,0.2,0,
                           0.3,0.2,1.5,0,
                           0,0,0,1), d, d),
                  matrix(c(1, -0.3, 0, 0,
                           -0.3, 1.2, 0.4, 0,
                           0, 0.4, 1, 0,
                           0, 0, 0, 0.8), d, d))

true_probs <- c(0.3, 0.4, 0.3)  # proportions

# Attribution des points aux clusters
z_true <- sample(1:K, size = n, replace = TRUE, prob = true_probs)

# Génération des points
X <- matrix(NA, nrow = n, ncol = d)
for(i in 1:n) {
  X[i,] <- mvrnorm(1, mu = true_means[[z_true[i]]], Sigma = true_covs[[z_true[i]]])
}

# ============================================================
# 2. Définition du modèle NIMBLE
# ============================================================
# On utilise un modèle de mélange gaussien :
# - pi ~ Dirichlet (proportions des clusters)
# - mu_k ~ Normale multidimensionnelle
# - Prec_k ~ Wishart (matrice de précision, inverse de la covariance)
# - z_i ~ Categorical(pi)
# - X_i ~ Normale multidimensionnelle selon le cluster z_i
# ============================================================

code <- nimbleCode({
  # Prior sur les proportions des clusters
  pi[1:K] ~ ddirch(alpha[1:K])
  
  # Priors sur les paramètres des clusters
  for(k in 1:K) {
    mu[k,1:d] ~ dmnorm(mean = mu0[1:d], cov = cov0[1:d,1:d]) # prior sur les moyennes
    Prec[k,1:d,1:d] ~ dwish(R[1:d,1:d], df)                  # précision (inverse covariance)
  }
  
  # Attribution des points et vraisemblance
  for(i in 1:N) {
    z[i] ~ dcat(pi[1:K])                                    # attribution au cluster
    X[i,1:d] ~ dmnorm(mean = mu[z[i],1:d], prec = Prec[z[i],1:d,1:d])
  }
})

# ============================================================
# 3. Constantes et données pour NIMBLE
# ============================================================

constants <- list(
  K = K,
  N = n,
  d = d,
  mu0 = rep(0, d),          # prior vague sur les moyennes
  cov0 = diag(d) * 10,      # variance large
  R = diag(d),              # matrice identité pour Wishart
  df = d + 1,               # degrés de liberté > d
  alpha = rep(1, K)         # prior uniforme sur les proportions
)

data <- list(X = X)

# Initialisation (z aléatoire)
inits <- list(z = sample(1:K, n, replace = TRUE))

# ============================================================
# 4. Compilation du modèle et MCMC
# ============================================================

model <- nimbleModel(code, data = data, constants = constants, inits = inits)
Cmodel <- compileNimble(model)

# Configuration du MCMC
conf <- configureMCMC(model, monitors = c("pi", "mu", "Prec", "z"))
Rmcmc <- buildMCMC(conf)
Cmcmc <- compileNimble(Rmcmc, project = model)

# Lancement du MCMC
samples <- runMCMC(Cmcmc, niter = 6000, nburnin = 1000, thin = 5, nchains = 1, setSeed = TRUE)

# ============================================================
# 5. Résultats : estimation des paramètres et attribution
# ============================================================

# Résumé des proportions de clusters
pi_idx <- grep("pi", colnames(samples))
pi_est <- colMeans(samples[, pi_idx])
cat("Proportions estimées (pi):", pi_est, "\n\n")

# Moyennes et covariances estimées
for (k in 1:K) {
  mu_idx <- grep(paste0("mu\\[",k,","), colnames(samples))
  mu_est <- colMeans(samples[, mu_idx])
  
  Prec_idx <- grep(paste0("Prec\\[",k,","), colnames(samples))
  Prec_est <- matrix(colMeans(samples[, Prec_idx]), d, d)
  Sigma_est <- solve(Prec_est)  # covariance estimée
  
  cat("=== Cluster", k, "===\n")
  cat("Moyenne estimée:\n")
  print(mu_est)
  cat("Covariance estimée:\n")
  print(Sigma_est)
  cat("\n")
}

# Attribution finale des points (MAP sur les z échantillonnés)
z_idx <- grep("z\\[", colnames(samples))
z_samples <- samples[, z_idx]

# Pour chaque point, attribution au cluster le plus fréquent dans l'échantillon
z_est <- apply(z_samples, 2, function(x) {
  as.numeric(names(sort(table(x), decreasing = TRUE)[1]))
})

# Afficher les 10 premières attributions
cat("Attribution des 10 premiers points :", z_est[1:10], "\n")
