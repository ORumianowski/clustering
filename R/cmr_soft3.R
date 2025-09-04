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
print(samples$summary)

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

