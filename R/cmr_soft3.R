library(nimble)
library(MASS)

# ----- PARAMÈTRES -----
set.seed(123)
N <- 150
K <- 3
d <- 4
Tmax <- 5

# Proportions et clusters
pi_true <- c(0.5, 0.3, 0.2)
z_true <- sample(1:K, N, replace=TRUE, prob=pi_true)

# Moyennes bien séparées
mu_true <- list(
  c(0, 0, 0, 0),
  c(4, 4, 4, 4), 
  c(-4, 4, -4, 4)
)

# Matrices de covariance DIAGONALES (évite la singularité)
Sigma_true <- list(
  diag(c(1.0, 1.0, 1.0, 1.0)),
  diag(c(1.5, 1.5, 1.5, 1.5)),
  diag(c(1.0, 2.0, 1.0, 2.0))
)

# Survie et détection
phi_true <- c(0.9, 0.7, 0.5)
p_true <- 0.6

# Coordonnées
X <- matrix(NA, nrow=N, ncol=d)
for(i in 1:N){
  k <- z_true[i]
  X[i,] <- mvrnorm(1, mu_true[[k]], Sigma_true[[k]])
}

# États latents et observations
state <- matrix(NA, nrow=N, ncol=Tmax)
y <- matrix(NA, nrow=N, ncol=Tmax)

for(i in 1:N){
  k <- z_true[i]
  state[i,1] <- 1
  for(t in 2:Tmax){
    if(state[i,t-1] == 1) {
      state[i,t] <- sample(1:2, size=1, prob=c(phi_true[k], 1 - phi_true[k]))
    } else {
      state[i,t] <- 2
    }
  }
  for(t in 1:Tmax){
    if(state[i,t] == 1){
      y[i,t] <- sample(1:2, size=1, prob=c(p_true, 1 - p_true))
    } else {
      y[i,t] <- 2
    }
  }
}

# ----- CODE NIMBLE SIMPLIFIÉ -----

code <- nimbleCode({
  # Priors
  pi[1:K] ~ ddirch(alpha[1:K])
  p ~ dbeta(1, 1)
  
  for(k in 1:K){
    # Moyennes avec prior plus informatif
    for(j in 1:d) {
      mu[k, j] ~ dnorm(0, var = 10)
    }
    
    # Précisons diagonales seulement (évite matrices singulières)
    for(j in 1:d) {
      sigma[k, j] ~ dinvgamma(shape = 2.1, scale = 1.1)  # shape > 2 pour éviter singularité
      prec[k, j] <- 1 / sigma[k, j]
    }
    
    phi[k] ~ dbeta(1, 1)
    
    # Matrice de transition
    transition[k, 1, 1] <- phi[k]
    transition[k, 1, 2] <- 1 - phi[k]
    transition[k, 2, 1] <- 0
    transition[k, 2, 2] <- 1
  }
  
  # Matrice d'observation
  obs[1, 1] <- p
  obs[1, 2] <- 1 - p
  obs[2, 1] <- 0
  obs[2, 2] <- 1
  
  for(i in 1:N){
    z[i] ~ dcat(pi[1:K])
    
    # Modèle pour X avec covariance diagonale
    for(j in 1:d) {
      X[i, j] ~ dnorm(mu[z[i], j], prec[z[i], j])
    }
    
    state[i, 1] ~ dcat(init[1:2])
    
    for(t in 2:Tmax){
      state[i, t] ~ dcat(transition[z[i], state[i, t-1], 1:2])
    }
    
    for(t in 1:Tmax){
      y[i, t] ~ dcat(obs[state[i, t], 1:2])
    }
  }
})

# ----- CONSTANTES -----

constants <- list(
  N = N,
  K = K,
  d = d,
  Tmax = Tmax,
  alpha = rep(1, K),
  init = c(1, 0)
)

data <- list(
  X = X,
  y = y
)

# ----- INITIALISATION -----

inits <- list(
  z = sample(1:K, N, replace = TRUE),
  mu = matrix(c(0, 0, 0, 0,
                3, 3, 3, 3,
                -3, 3, -3, 3), nrow = K, ncol = d, byrow = TRUE),
  sigma = matrix(rep(1, K*d), nrow = K, ncol = d),
  phi = c(0.8, 0.6, 0.4),
  p = 0.5,
  state = matrix(1, N, Tmax)
)

# ----- VÉRIFICATION -----

model <- nimbleModel(code, constants = constants, data = data, inits = inits)
cat("Log-probabilité initiale:", model$calculate(), "\n")

# ----- MCMC -----

Cmodel <- compileNimble(model)

conf <- configureMCMC(model, monitors = c("pi", "mu", "sigma", "phi", "p"))
Rmcmc <- buildMCMC(conf)
Cmcmc <- compileNimble(Rmcmc, project = model)

# Exécution
samples <- runMCMC(Cmcmc, 
                   niter = 3000, 
                   nburnin = 1000, 
                   thin = 5,
                   progressBar = TRUE)

# ----- RÉSULTATS -----
# print(samples$summary)

library(ggplot2)
library(tidyr)
library(dplyr)

# ---- extraction des échantillons bruts ----
samples_df <- as.data.frame(samples)  

# ---- valeurs vraies pour TOUS les paramètres ----
true_values <- c(
  # mu
  "mu[1, 1]" = 0,  "mu[1, 2]" = 0,  "mu[1, 3]" = 0,  "mu[1, 4]" = 0,
  "mu[2, 1]" = 4,  "mu[2, 2]" = 4,  "mu[2, 3]" = 4,  "mu[2, 4]" = 4,
  "mu[3, 1]" = -4, "mu[3, 2]" = 4,  "mu[3, 3]" = -4, "mu[3, 4]" = 4,
  
  # p
  "p" = 0.6,
  
  # phi
  "phi[1]" = 0.9, "phi[2]" = 0.7, "phi[3]" = 0.5,
  
  # pi
  "pi[1]" = 0.3, "pi[2]" = 0.4, "pi[3]" = 0.3,
  
  # sigma (diagonales des covariances)
  "sigma[1, 1]" = 1.0, "sigma[1, 2]" = 1.0, "sigma[1, 3]" = 1.0, "sigma[1, 4]" = 1.0,
  "sigma[2, 1]" = 1.5, "sigma[2, 2]" = 1.5, "sigma[2, 3]" = 1.5, "sigma[2, 4]" = 1.5,
  "sigma[3, 1]" = 1.0, "sigma[3, 2]" = 2.0, "sigma[3, 3]" = 1.0, "sigma[3, 4]" = 2.0
)

# ---- passage en format long ----
df_long <- samples_df %>%
  pivot_longer(cols = everything(),
               names_to = "param", values_to = "value") %>%
  mutate(true_value = true_values[param],
         group = case_when(
           grepl("^mu", param) ~ "mu",
           grepl("^sigma", param) ~ "sigma",
           grepl("^phi", param) ~ "phi",
           grepl("^pi", param) ~ "pi",
           param == "p" ~ "p"
         ))

# ---- fonction générique de plot ----
plot_params <- function(df, group_name) {
  df %>%
    filter(group == group_name) %>%
    ggplot(aes(x = value)) +
    geom_histogram(aes(y = ..density..), bins = 30,
                   fill = "skyblue", color = "white", alpha = 0.7) +
    geom_vline(aes(xintercept = true_value),
               color = "red", linetype = "dashed", size = 1) +
    facet_wrap(~ param, scales = "free", ncol = 4) +
    theme_minimal(base_size = 14) +
    labs(title = paste("Distributions postérieures des paramètres -", group_name),
         x = "Valeur", y = "Densité")
}

# ---- plots par type ----
plot_params(df_long, "mu")
plot_params(df_long, "sigma")
plot_params(df_long, "phi")
plot_params(df_long, "pi")
plot_params(df_long, "p")


# ---- ALGORITHME DE RÉALIGNEMENT DES CLUSTERS ----

# Fonction pour trouver la meilleure permutation des labels
find_best_permutation <- function(mu_samples, mu_true) {
  K <- nrow(mu_true)
  permutations <- combinat::permn(1:K)  # toutes les permutations possibles
  
  best_perm <- NULL
  best_error <- Inf
  
  for(perm in permutations) {
    error <- sum((mu_samples[perm,] - mu_true)^2)
    if(error < best_error) {
      best_error <- error
      best_perm <- perm
    }
  }
  
  return(best_perm)
}

# Calcul des moyennes médianes par cluster
mu_est <- matrix(NA, nrow = K, ncol = d)
for(k in 1:K) {
  for(j in 1:d) {
    mu_est[k, j] <- median(samples_df[, paste0("mu[", k, ", ", j, "]")])
  }
}

# Trouver la meilleure permutation
best_permutation <- find_best_permutation(mu_est, do.call(rbind, mu_true))
cat("Meilleure permutation des labels:", best_permutation, "\n")

# ---- RÉALIGNEMENT DES VALEURS VRAIES ----

# Réorganiser les vraies valeurs selon la permutation
true_values_aligned <- c(
  # mu
  "mu[1, 1]" = mu_true[[1]][1], 
  "mu[1, 2]" = mu_true[[1]][2],
  "mu[1, 3]" = mu_true[[1]][3],
  "mu[1, 4]" = mu_true[[1]][4],
  "mu[2, 1]" = mu_true[[2]][1],
  "mu[2, 2]" = mu_true[[2]][2],
  "mu[2, 3]" = mu_true[[2]][3],
  "mu[2, 4]" = mu_true[[2]][4],
  "mu[3, 1]" = mu_true[[3]][1],
  "mu[3, 2]" = mu_true[[3]][2],
  "mu[3, 3]" = mu_true[[3]][3],
  "mu[3, 4]" = mu_true[[3]][4],
  
  # p
  "p" = 0.6,
  
  # phi
  "phi[1]" = phi_true[1],
  "phi[2]" = phi_true[2],
  "phi[3]" = phi_true[3],
  
  # pi
  "pi[1]" = pi_true[1],
  "pi[2]" = pi_true[2],
  "pi[3]" = pi_true[3],
  
  # sigma
  "sigma[1, 1]" = Sigma_true[[1]][1,1],
  "sigma[1, 2]" = Sigma_true[[1]][2,2],
  "sigma[1, 3]" = Sigma_true[[1]][3,3],
  "sigma[1, 4]" = Sigma_true[[1]][4,4],
  "sigma[2, 1]" = Sigma_true[[2]][1,1],
  "sigma[2, 2]" = Sigma_true[[2]][2,2],
  "sigma[2, 3]" = Sigma_true[[2]][3,3],
  "sigma[2, 4]" = Sigma_true[[2]][4,4],
  "sigma[3, 1]" = Sigma_true[[3]][1,1],
  "sigma[3, 2]" = Sigma_true[[3]][2,2],
  "sigma[3, 3]" = Sigma_true[[3]][3,3],
  "sigma[3, 4]" = Sigma_true[[3]][4,4]
)


# ---- RÉALIGNEMENT DES ÉCHANTILLONS MCMC ----

# Créer une nouvelle dataframe avec les échantillons réalignés
samples_aligned <- samples_df

# Réaligner les paramètres de cluster
for(k in 1:K) {
  new_k <- which(best_permutation == k)  # quel cluster devient k
  
  # Réaligner mu
  for(j in 1:d) {
    old_name <- paste0("mu[", new_k, ", ", j, "]")
    new_name <- paste0("mu[", k, ", ", j, "]")
    samples_aligned[[new_name]] <- samples_df[[old_name]]
  }
  
  # Réaligner phi
  old_phi <- paste0("phi[", new_k, "]")
  new_phi <- paste0("phi[", k, "]")
  samples_aligned[[new_phi]] <- samples_df[[old_phi]]
  
  # Réaligner pi
  old_pi <- paste0("pi[", new_k, "]")
  new_pi <- paste0("pi[", k, "]")
  samples_aligned[[new_pi]] <- samples_df[[old_pi]]
  
  # Réaligner sigma
  for(j in 1:d) {
    old_sigma <- paste0("sigma[", new_k, ", ", j, "]")
    new_sigma <- paste0("sigma[", k, ", ", j, "]")
    samples_aligned[[new_sigma]] <- samples_df[[old_sigma]]
  }
}

# ---- VISUALISATION AVEC CLUSTERS RÉALIGNÉS ----

df_long_aligned <- samples_aligned %>%
  pivot_longer(cols = everything(),
               names_to = "param", values_to = "value") %>%
  mutate(true_value = true_values_aligned[param],
         group = case_when(
           grepl("^mu", param) ~ "mu",
           grepl("^sigma", param) ~ "sigma",
           grepl("^phi", param) ~ "phi",
           grepl("^pi", param) ~ "pi",
           param == "p" ~ "p"
         ))

# ---- plots par type avec clusters réalignés ----
plot_params(df_long_aligned, "mu") +
  ggtitle("Distributions postérieures des moyennes (clusters réalignés)")

plot_params(df_long_aligned, "sigma") +
  ggtitle("Distributions postérieures des variances (clusters réalignés)")

plot_params(df_long_aligned, "phi") +
  ggtitle("Distributions postérieures des survies (clusters réalignés)")

plot_params(df_long_aligned, "pi") +
  ggtitle("Distributions postérieures des proportions (clusters réalignés)")

plot_params(df_long_aligned, "p") +
  ggtitle("Distribution postérieure de la probabilité de détection")

# ---- ÉVALUATION QUANTITATIVE ----
cat("=== PERFORMANCE DU MODÈLE ===\n")

# Calcul des biais
estimates <- df_long_aligned %>%
  group_by(param) %>%
  summarise(median_est = median(value),
            true_val = first(true_value)) %>%
  mutate(bias = median_est - true_val,
         relative_bias = abs(bias)/true_val * 100)

print(estimates)

cat("\Biais moyen absolu:", mean(abs(estimates$bias)), "\n")
cat("Biais relatif moyen:", mean(estimates$relative_bias, na.rm = TRUE), "%\n")

