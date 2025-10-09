
# Script complet : simulation, 3 modèles (integrated, survival-only ordered phi with fixed pi, connectivity-only),


library(nimble)
library(MASS)
library(coda)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(cowplot)
library(gtools)

set.seed(123)

# --------------------------
# 1) SIMULATION (identique)
# --------------------------
# Nombre d'individus total
N <- 150
# Nombre de clusters
K <- 3
# Nombre de dimensions au cluster 
d <- 4 # (2 coordonées en breeding, 2 coordonnées en wintering)
# Nombre d'occasion de capture
Tmax <- 5

## CONNECTIVITE
# Repartion des individus echantilonnés dans les clusters
pi_true <- c(0.5, 0.3, 0.2)
# Assignation des individus à un cluster
cl_true <- sample(1:K, N, replace = TRUE, prob = pi_true)

# Coordonnnées des barycentres de clusters
# Les deux premiers sont proches proches mais pourront être différencer grâce à des survie très différentes
mu_true <- list(
  c(0, 0, 0, 0),
  c(1, 2, 0, 0),
  c(-4, 4, -4, 4)
)
# Covariance des gaussiennes
Sigma_true <- list(
  diag(c(1.0, 1.0, 1.0, 1.0)),
  diag(c(1.5, 1.5, 1.5, 1.5)),
  diag(c(1.0, 2.0, 1.0, 2.0))
)

# SURVIE
# Survie par cluster
phi_true <- c(0.8, 0.6, 0.5) 
# Probabilite de dectection (CMR)
p_true <- 0.6

# Generer les 4 coordonnées pour tous les individus sachant son cluster
X <- matrix(NA, nrow = N, ncol = d)
for (i in 1:N) {
  k <- cl_true[i]
  X[i, ] <- mvrnorm(1, mu_true[[k]], Sigma_true[[k]])
}

# Etats & observations y
state <- matrix(NA, nrow = N, ncol = Tmax)
y <- matrix(NA, nrow = N, ncol = Tmax)

# Single batch: tous les individus sont capturés à la premiere occasion
for (i in 1:N) {
  k <- cl_true[i]
  state[i,1] <- 1
  # States
  for (t in 2:Tmax) {
    if (state[i,t-1] == 1) state[i,t] <- sample(1:2, 1, prob = c(phi_true[k], 1-phi_true[k])) else state[i,t] <- 2
  }
  # Obs
  for (t in 1:Tmax) {
    if (state[i,t] == 1) y[i,t] <- sample(1:2, 1, prob = c(p_true, 1-p_true)) else y[i,t] <- 2
  }
}

# --------------------------
# 2) CODES NIMBLE
# --------------------------

# 2.a) Intégré
code_integrated <- nimbleCode({
  # Probabilite d appartenir a un cluster
  pi[1:K] ~ ddirch(alpha[1:K])
  # Paramètres des gaussiènnes
  for(k in 1:K){
    for(j in 1:d){
      # Barycentres des gaussiennes
      mu[k,j] ~ dnorm(0, var = 10)
      # Covariance des gaussiennes
      sigma[k,j] ~ dinvgamma(shape = 2.1, scale = 1.1)
      # Il est conseillé de travailler avec la precision, plutot que le covariance
      prec[k,j] <- 1 / sigma[k,j]
    }
    # Survie
    phi[k] ~ dbeta(1,1)
    # Probabilite de detection
    p ~ dbeta(1,1)
    # Matrice de transition d'état
    transition[k,1,1] <- phi[k]; transition[k,1,2] <- 1-phi[k]
    transition[k,2,1] <- 0; transition[k,2,2] <- 1
  }
  # Matrice de transition d'état latent vers état observé
  obs[1,1] <- p; obs[1,2] <- 1-p; obs[2,1] <- 0; obs[2,2] <- 1
  
  # Processus de transition
  for(i in 1:N){
    # Vraisemblance de la connectivité
    # Assignation d'un cluster à l'individu
    cl[i] ~ dcat(pi[1:K])
    # Sachant son cluster, détermination de la probabilité d'avoir cette position sacant la distribution du cluster
    for(j in 1:d) X[i,j] ~ dnorm(mu[cl[i],j], prec[cl[i],j])
    
    # Vraisemblance de la survie
    state[i,1] <- 1 #  ~ dcat(init[1:2])
    for(t in 2:Tmax) state[i,t] ~ dcat(transition[cl[i], state[i,t-1], 1:2])
    for(t in 1:Tmax) y[i,t] ~ dcat(obs[state[i,t], 1:2])
  }
})

# 2.b) Survie seule : pi fixé, phi ordonné via delta
code_survival_ordered <- nimbleCode({
  # pi_fix fourni en data (pi_true)
  p ~ dbeta(1,1)
  # les delta sont des termes accesoires pour ordoner les survies
  for(k in 1:K) delta[k] ~ dbeta(1,1)
  # construire phi croissants
  phi[1] <- delta[1]
  for(k in 2:K) phi[k] <- phi[k-1] + delta[k] * (1 - phi[k-1])
  for (k in 1:K) {
    transition[k,1,1] <- phi[k]; transition[k,1,2] <- 1-phi[k]
    transition[k,2,1] <- 0; transition[k,2,2] <- 1
  }
  obs[1,1] <- p; obs[1,2] <- 1-p; obs[2,1] <- 0; obs[2,2] <- 1
  for(i in 1:N){
    cl[i] ~ dcat(pi_fix[1:K])    # pi_fix fourni comme data = pi_true
    state[i,1] ~ dcat(init[1:2])
    for(t in 2:Tmax) state[i,t] ~ dcat(transition[cl[i], state[i,t-1], 1:2])
    for(t in 1:Tmax) y[i,t] ~ dcat(obs[state[i,t], 1:2])
  }
})

# 2.c) Connectivité seule
code_connectivity_only <- nimbleCode({
  pi[1:K] ~ ddirch(alpha[1:K])
  for(k in 1:K){
    for(j in 1:d){
      mu[k,j] ~ dnorm(0, var=10)
      sigma[k,j] ~ dinvgamma(shape=2.1, scale=1.1)
      prec[k,j] <- 1/sigma[k,j]
    }
  }
  for(i in 1:N){
    cl[i] ~ dcat(pi[1:K])
    for(j in 1:d) X[i,j] ~ dnorm(mu[cl[i],j], prec[cl[i],j])
  }
})

# --------------------------
# 3) CONSTANTES, DATA, INITS
# --------------------------
constants <- list(N=N, K=K, d=d, Tmax=Tmax, alpha=rep(1,K), init=c(1,0))
data_integrated <- list(X = X, y = y)
data_conn_only <- list(X = X)
data_surv_only <- list(y = y,
                       pi_fix = pi_true) # on aide le modèle de survie seul en lui donnant directement la proportion d'individus par cluster.


inits_common <- list(
  cl = sample(1:K, N, replace = TRUE),   
  mu = matrix(c(0,0,0,0, 3,3,3,3, -3,3,-3,3), nrow = K, byrow = TRUE),
  sigma = matrix(1, nrow = K, ncol = d),
  delta = rep(0.5, K),
  p = 0.5,
  state = matrix(1, nrow = N, ncol = Tmax)
)

inits_integrated <- inits_common

inits_conn <- list(
  cl = inits_common$cl,   
  mu = inits_common$mu,
  sigma = inits_common$sigma,
  pi = rep(1/K, K)
)

inits_surv <- list(
  cl = sample(1:K, N, replace = TRUE),  
  delta = rep(0.6, K),
  p = 0.5,
  state = matrix(1, N, Tmax)
)


# --------------------------
# 4) RUN NIMBLE helper
# --------------------------
run_nimble <- function(code, constants, data, inits, monitors,
                       niter = 6000, nburnin = 2000, thin = 4, progressBar=TRUE) {
  model <- nimbleModel(code, constants = constants, data = data, inits = inits, calculate = TRUE)
  cmodel <- compileNimble(model)
  conf <- configureMCMC(model, monitors = monitors)
  mcmc <- buildMCMC(conf)
  cmcmc <- compileNimble(mcmc, project = model)
  samp <- runMCMC(cmcmc, niter = niter, nburnin = nburnin, thin = thin, progressBar = progressBar)
  return(list(model = model, cmodel = cmodel, samples = as.mcmc(samp)))
}

# --------------------------
# 5) LANCER LES CHAÎNES (ajuster niter si nécessaire)
# --------------------------
niter <- 6000; nburnin <- 2000; thin <- 4

cat("Running integrated model...\n")
res_int <- run_nimble(code_integrated, constants, data_integrated, inits_integrated,
                      monitors = c("pi","mu","sigma","phi","p"),
                      niter=niter, nburnin=nburnin, thin=thin)

cat("Running survival-only (pi fixed, phi ordered)...\n")
res_surv <- run_nimble(code_survival_ordered, constants, data_surv_only, inits_surv,
                       monitors = c("phi","delta","p"),
                       niter=niter, nburnin=nburnin, thin=thin)

cat("Running connectivity-only model...\n")
res_conn <- run_nimble(code_connectivity_only, constants, data_conn_only, inits_conn,
                       monitors = c("pi","mu","sigma"),
                       niter=niter, nburnin=nburnin, thin=thin)

# --------------------------
# 6) Convertir samples en data.frame
# --------------------------
samples_int_df <- as.data.frame(as.matrix(res_int$samples))
samples_surv_df <- as.data.frame(as.matrix(res_surv$samples))
samples_conn_df <- as.data.frame(as.matrix(res_conn$samples))

# --------------------------
# 7) Fonctions d'extraction & permutation
# --------------------------

# Les clusters n'ont pas nécessairement la même étiquette entre elle et avec les valeurs simulées. 
# Il faut identifier la combianaisons la plus probable.

# Extraire matrice K x d de médianes pour mu 
extract_mu_median_mat <- function(samples_df) {
  mu_cols <- grep("^mu\\[", names(samples_df), value = TRUE)
  if(length(mu_cols) == 0) return(NULL)
  # parse names mu[k,j] (nimble prints "mu[1,1]" or "mu[1, 1]" depending)
  mu_mat <- matrix(NA, nrow = K, ncol = d)
  for(col in mu_cols) {
    nums <- as.integer(regmatches(col, gregexpr("\\d+", col))[[1]])
    k <- nums[1]; j <- nums[2]
    mu_mat[k, j] <- median(samples_df[[col]])
  }
  return(mu_mat)
}

# Extraire vecteur phi median
extract_phi_median_vec <- function(samples_df) {
  phi_cols <- grep("^phi\\[", names(samples_df), value = TRUE)
  if(length(phi_cols) == 0) return(NULL)
  phi_vec <- numeric(K)
  for(col in phi_cols) {
    k <- as.integer(regmatches(col, gregexpr("\\d+", col))[[1]])
    phi_vec[k] <- median(samples_df[[col]])
  }
  return(phi_vec)
}

# Trouver la meilleure permutation pour matrices (mu)
find_best_perm_matrix <- function(est_mat, true_mat) {
  perms <- permutations(n = K, r = K, v = 1:K)
  best_perm <- NULL; best_err <- Inf
  for(i in seq_len(nrow(perms))) {
    perm <- perms[i, ]
    # we adopt convention: perm[new_k] = old_k  (est_mat[old_k,] -> true row new_k)
    err <- sum((est_mat[perm, , drop = FALSE] - true_mat)^2)
    if(err < best_err) { best_err <- err; best_perm <- perm }
  }
  return(best_perm)
}

# Trouver permutation pour vector (phi)
find_best_perm_vector <- function(est_vec, true_vec) {
  perms <- permutations(n = K, r = K, v = 1:K)
  best_perm <- NULL; best_err <- Inf
  for(i in seq_len(nrow(perms))) {
    perm <- perms[i, ]
    err <- sum((est_vec[perm] - true_vec)^2)
    if(err < best_err) { best_err <- err; best_perm <- perm }
  }
  return(best_perm)
}

# Fonction pour renommer les colonnes d'un dataframe de samples selon permutation perm
# perm est tel que perm[new_k] = old_k  (c.f. above)
rename_samples_by_perm <- function(samples_df, perm) {
  inv_perm <- integer(K)  # inv_perm[old_k] = new_k
  for(newk in 1:K) inv_perm[perm[newk]] <- newk
  new_names <- names(samples_df)
  for(i in seq_along(new_names)) {
    nm <- new_names[i]
    # mu
    if(grepl("^mu\\[", nm)) {
      nums <- as.integer(regmatches(nm, gregexpr("\\d+", nm))[[1]])
      oldk <- nums[1]; j <- nums[2]
      newk <- inv_perm[oldk]
      new_names[i] <- paste0("mu[", newk, ",", j, "]")
      next
    }
    # sigma
    if(grepl("^sigma\\[", nm)) {
      nums <- as.integer(regmatches(nm, gregexpr("\\d+", nm))[[1]])
      oldk <- nums[1]; j <- nums[2]
      newk <- inv_perm[oldk]
      new_names[i] <- paste0("sigma[", newk, ",", j, "]")
      next
    }
    # phi
    if(grepl("^phi\\[", nm)) {
      oldk <- as.integer(regmatches(nm, gregexpr("\\d+", nm))[[1]])
      newk <- inv_perm[oldk]
      new_names[i] <- paste0("phi[", newk, "]")
      next
    }
    # pi
    if(grepl("^pi\\[", nm)) {
      oldk <- as.integer(regmatches(nm, gregexpr("\\d+", nm))[[1]])
      newk <- inv_perm[oldk]
      new_names[i] <- paste0("pi[", newk, "]")
      next
    }
    # else keep same
  }
  names(samples_df) <- new_names
  return(samples_df)
}

# --------------------------
# 8) Calculer permutations & réordonner samples pour chaque modèle
# --------------------------
mu_true_mat <- do.call(rbind, mu_true)  # K x d

# integrated: find perm using mu
mu_int_est <- extract_mu_median_mat(samples_int_df)
if(!is.null(mu_int_est)) {
  perm_int <- find_best_perm_matrix(mu_int_est, mu_true_mat)
  cat("Permutation optimale (integrated):", perm_int, "\n")
  samples_int_df_aligned <- rename_samples_by_perm(samples_int_df, perm_int)
} else {
  samples_int_df_aligned <- samples_int_df
  perm_int <- 1:K
}

# connectivity-only: find perm using mu
mu_conn_est <- extract_mu_median_mat(samples_conn_df)
if(!is.null(mu_conn_est)) {
  perm_conn <- find_best_perm_matrix(mu_conn_est, mu_true_mat)
  cat("Permutation optimale (connectivity-only):", perm_conn, "\n")
  samples_conn_df_aligned <- rename_samples_by_perm(samples_conn_df, perm_conn)
} else {
  samples_conn_df_aligned <- samples_conn_df
  perm_conn <- 1:K
}

# survival-only: align using phi (we have ordered phi inside model but order refers to model phi[1]<phi[2]<phi[3];
# we want to map these ordered phi to the true phi vector)
phi_surv_est <- extract_phi_median_vec(samples_surv_df)
perm_surv <- find_best_perm_vector(phi_surv_est, phi_true)
cat("Permutation optimale (survival-only):", perm_surv, "\n")
samples_surv_df_aligned <- rename_samples_by_perm(samples_surv_df, perm_surv)

# --------------------------
# 9) Construire long dataframe des paramètres alignés
# --------------------------
df_long_from_samples <- function(samples_df, method_label) {
  df <- as.data.frame(samples_df) %>% mutate(.iter = row_number(), method = method_label) %>%
    pivot_longer(cols = -c(.iter, method), names_to = "param", values_to = "value")
  df <- df %>%
    mutate(type = case_when(
      grepl("^mu\\[", param) ~ "mu",
      grepl("^sigma\\[", param) ~ "sigma",
      grepl("^phi\\[", param) ~ "phi",
      grepl("^pi\\[", param) ~ "pi",
      param == "p" ~ "p",
      TRUE ~ "other"
    )) %>%
    rowwise() %>%
    mutate(nums = list(as.integer(regmatches(param, gregexpr("\\d+", param))[[1]])),
           cluster = ifelse(length(nums) >= 1, nums[1], NA_integer_),
           dim = ifelse(length(nums) >= 2, nums[2], NA_integer_)) %>%
    select(-nums)
  return(df)
}

df_int_long <- df_long_from_samples(samples_int_df_aligned, "integrated")
df_conn_long <- df_long_from_samples(samples_conn_df_aligned, "connectivity_only")
df_surv_long <- df_long_from_samples(samples_surv_df_aligned, "survival_only")

all_long <- bind_rows(df_int_long, df_conn_long, df_surv_long)

# --------------------------
# 10) RÉSUMÉS POSTERIEURS (moyenne, IC95%) et joindre vraies valeurs
# --------------------------
summary_by_param <- all_long %>%
  filter(type %in% c("mu","phi","sigma","pi","p")) %>%
  group_by(method, type, param, cluster, dim) %>%
  summarise(mean = mean(value), sd = sd(value),
            lower = quantile(value, .025), upper = quantile(value, .975),
            width = upper - lower, .groups = "drop")

# table des vraies valeurs
true_df_mu <- expand.grid(cluster = 1:K, dim = 1:d) %>%
  mutate(true = map2_dbl(cluster, dim, ~ mu_true[[.x]][.y]), type = "mu")
true_df_phi <- data.frame(cluster = 1:K, true = phi_true, type = "phi")
true_df_sigma <- expand.grid(cluster = 1:K, dim = 1:d) %>%
  mutate(true = map2_dbl(cluster, dim, ~ Sigma_true[[.x]][.y,.y]), type = "sigma")
true_df_pi <- data.frame(cluster = 1:K, true = pi_true, type = "pi")
true_df <- bind_rows(true_df_mu, true_df_phi, true_df_sigma, true_df_pi)

summary_with_true <- summary_by_param %>%
  left_join(true_df, by = c("cluster","type")) %>%
  mutate(bias = mean - true, rmse = sqrt((mean-true)^2 + sd^2))

# --------------------------
# 11) GRAPHIQUES (après alignement)
# --------------------------
method_colors <- c("integrated" = "#1b9e77", "survival_only" = "#d95f02", "connectivity_only" = "#7570b3")

# density phi
p_phi_density <- all_long %>% filter(type == "phi") %>%
  ggplot(aes(x = value, color = method, fill = method)) +
  geom_density(alpha = .25) +
  facet_wrap(~ cluster, scales = "free") +
  geom_vline(data = data.frame(cluster = 1:K, true = phi_true), aes(xintercept = true), linetype = "dashed") +
  scale_color_manual(values = method_colors) + scale_fill_manual(values = method_colors) +
  labs(title = "Densités postérieures de phi", x = expression(phi), y = "Densité") +
  theme_minimal()

# density mu
p_mu_density <- all_long %>% filter(type == "mu") %>%
  ggplot(aes(x = value, color = method, fill = method)) +
  geom_density(alpha = .25) +
  facet_grid(cluster ~ dim, scales = "free") +
  geom_vline(data = true_df_mu, aes(xintercept = true), linetype = "dashed") +
  scale_color_manual(values = method_colors) + scale_fill_manual(values = method_colors) +
  labs(title = "Densités postérieures de mu", x = "mu", y = "Densité") +
  theme_minimal()

# --------------------------
# 12) Plot
# --------------------------

print(p_phi_density); print(p_mu_density)

