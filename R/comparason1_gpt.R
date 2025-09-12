# analysis_integrated_vs_separate.R
# Script complet : simulation (déjà faite), 3 modèles (intégré, survie seule, connectivité seule),
# extraction MCMC, diagnostics, graphiques comparatifs.
#
# Dépendances : nimble, MASS, coda, ggplot2, dplyr, tidyr, purrr, cowplot
# Installer si nécessaire :
# install.packages(c("MASS","coda","ggplot2","dplyr","tidyr","purrr","cowplot"))
# install.packages("nimble") # si pas déjà installé

library(nimble)
library(MASS)
library(coda)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(cowplot)

set.seed(123)

# --------------------------
# 1) SIMULATION (ton code)
# --------------------------
N <- 150
K <- 3
d <- 4
Tmax <- 5

pi_true <- c(0.5, 0.3, 0.2)
z_true <- sample(1:K, N, replace = TRUE, prob = pi_true)

mu_true <- list(
  c(0, 0, 0, 0),
  c(4, 4, 4, 4),
  c(-4, 4, -4, 4)
)

Sigma_true <- list(
  diag(c(1.0, 1.0, 1.0, 1.0)),
  diag(c(1.5, 1.5, 1.5, 1.5)),
  diag(c(1.0, 2.0, 1.0, 2.0))
)

phi_true <- c(0.9, 0.7, 0.5)
p_true <- 0.6

# Coordonnées X (4 colonnes : breeding1, breeding2, wintering1, wintering2 per ta spécification)
X <- matrix(NA, nrow = N, ncol = d)
for (i in 1:N) {
  k <- z_true[i]
  X[i, ] <- mvrnorm(1, mu_true[[k]], Sigma_true[[k]])
}

# Etats latents et observations y
state <- matrix(NA, nrow = N, ncol = Tmax)
y <- matrix(NA, nrow = N, ncol = Tmax)

for (i in 1:N) {
  k <- z_true[i]
  state[i, 1] <- 1
  for (t in 2:Tmax) {
    if (state[i, t - 1] == 1) {
      state[i, t] <- sample(1:2, size = 1, prob = c(phi_true[k], 1 - phi_true[k]))
    } else {
      state[i, t] <- 2
    }
  }
  for (t in 1:Tmax) {
    if (state[i, t] == 1) {
      y[i, t] <- sample(1:2, size = 1, prob = c(p_true, 1 - p_true))
    } else {
      y[i, t] <- 2
    }
  }
}

# --------------------------
# 2) CODES NIMBLE POUR LES 3 MODÈLES
# --------------------------

# 2.a) Modèle intégré (ton code légèrement nettoyé)
code_integrated <- nimbleCode({
  # Priors
  pi[1:K] ~ ddirch(alpha[1:K])
  p ~ dbeta(1, 1)
  for (k in 1:K) {
    for (j in 1:d) {
      mu[k, j] ~ dnorm(0, var = 10)
      sigma[k, j] ~ dinvgamma(shape = 2.1, scale = 1.1)
      prec[k, j] <- 1 / sigma[k, j]
    }
    phi[k] ~ dbeta(1, 1)
    transition[k, 1, 1] <- phi[k]
    transition[k, 1, 2] <- 1 - phi[k]
    transition[k, 2, 1] <- 0
    transition[k, 2, 2] <- 1
  }
  obs[1, 1] <- p
  obs[1, 2] <- 1 - p
  obs[2, 1] <- 0
  obs[2, 2] <- 1
  
  for (i in 1:N) {
    z[i] ~ dcat(pi[1:K])
    for (j in 1:d) {
      X[i, j] ~ dnorm(mu[z[i], j], prec[z[i], j])
    }
    state[i, 1] ~ dcat(init[1:2])
    for (t in 2:Tmax) {
      state[i, t] ~ dcat(transition[z[i], state[i, t - 1], 1:2])
    }
    for (t in 1:Tmax) {
      y[i, t] ~ dcat(obs[state[i, t], 1:2])
    }
  }
})

# 2.b) Modèle survie seule
# On laisse pi comme paramètre à estimer (aide pour z), mais X n'est pas utilisé.
code_survival_only <- nimbleCode({
  pi[1:K] ~ ddirch(alpha[1:K])
  p ~ dbeta(1, 1)
  for (k in 1:K) {
    phi[k] ~ dbeta(1, 1)
    transition[k, 1, 1] <- phi[k]
    transition[k, 1, 2] <- 1 - phi[k]
    transition[k, 2, 1] <- 0
    transition[k, 2, 2] <- 1
  }
  obs[1, 1] <- p
  obs[1, 2] <- 1 - p
  obs[2, 1] <- 0
  obs[2, 2] <- 1
  
  for (i in 1:N) {
    z[i] ~ dcat(pi[1:K])
    state[i, 1] ~ dcat(init[1:2])
    for (t in 2:Tmax) {
      state[i, t] ~ dcat(transition[z[i], state[i, t - 1], 1:2])
    }
    for (t in 1:Tmax) {
      y[i, t] ~ dcat(obs[state[i, t], 1:2])
    }
  }
})

# 2.c) Modèle connectivité seule (X seulement)
code_connectivity_only <- nimbleCode({
  pi[1:K] ~ ddirch(alpha[1:K])
  for (k in 1:K) {
    for (j in 1:d) {
      mu[k, j] ~ dnorm(0, var = 10)
      sigma[k, j] ~ dinvgamma(shape = 2.1, scale = 1.1)
      prec[k, j] <- 1 / sigma[k, j]
    }
  }
  for (i in 1:N) {
    z[i] ~ dcat(pi[1:K])
    for (j in 1:d) {
      X[i, j] ~ dnorm(mu[z[i], j], prec[z[i], j])
    }
  }
})

# --------------------------
# 3) CONSTANTES, DATA, INITS
# --------------------------
constants <- list(N = N, K = K, d = d, Tmax = Tmax, alpha = rep(1, K), init = c(1, 0))
data_integrated <- list(X = X, y = y)
data_surv_only <- list(y = y)
data_conn_only <- list(X = X)

inits_common <- list(
  z = sample(1:K, N, replace = TRUE),
  mu = matrix(c(0, 0, 0, 0,
                3, 3, 3, 3,
                -3, 3, -3, 3), nrow = K, ncol = d, byrow = TRUE),
  sigma = matrix(rep(1, K * d), nrow = K, ncol = d),
  phi = c(0.8, 0.6, 0.4),
  p = 0.5,
  state = matrix(1, N, Tmax)
)

# Inits adaptés pour chaque modèle (nimble exige la liste inits de paramètres présents)
inits_integrated <- inits_common
inits_surv <- inits_common[c("z", "phi", "p", "state")]
inits_surv$pi <- rep(1 / K, K)
inits_conn <- inits_common[c("z", "mu", "sigma")]
inits_conn$pi <- rep(1 / K, K)

# --------------------------
# 4) FONCTION D'APPEL NIMBLE (compilation + run)
# --------------------------
run_nimble <- function(code, constants, data, inits, monitors,
                       niter = 3000, nburnin = 1000, thin = 2, showProgress = TRUE) {
  model <- nimbleModel(code, constants = constants, data = data, inits = inits, calculate = TRUE)
  cmodel <- compileNimble(model)
  conf <- configureMCMC(model, monitors = monitors)
  mcmc <- buildMCMC(conf)
  cmcmc <- compileNimble(mcmc, project = model)
  samples <- runMCMC(cmcmc, niter = niter, nburnin = nburnin, thin = thin, progressBar = showProgress)
  mcmc_coda <- as.mcmc(samples)
  return(list(model = model, cmodel = cmodel, mcmc = mcmc, samples = mcmc_coda))
}

# --------------------------
# 5) LANCEMENT DES CHAÎNES (adapter niter si lent)
# --------------------------
# Paramètres MCMC — augmente si tu veux plus de précision
niter <- 6000
nburnin <- 2000
thin <- 4

cat("Running integrated model...\n")
res_integrated <- run_nimble(code_integrated, constants, data_integrated, inits_integrated,
                             monitors = c("pi", "mu", "sigma", "phi", "p"),
                             niter = niter, nburnin = nburnin, thin = thin)

cat("Running survival-only model...\n")
res_surv <- run_nimble(code_survival_only, constants, data_surv_only, inits_surv,
                       monitors = c("pi", "phi", "p"),
                       niter = niter, nburnin = nburnin, thin = thin)

cat("Running connectivity-only model...\n")
res_conn <- run_nimble(code_connectivity_only, constants, data_conn_only, inits_conn,
                       monitors = c("pi", "mu", "sigma"),
                       niter = niter, nburnin = nburnin, thin = thin)

# --------------------------
# 6) EXTRACTION ET TRANSFORMATION DES ÉCHANTILLONS
# --------------------------

# Helper pour convertir as.data.frame des mcmc
mcmc_to_df <- function(mcmc_obj) {
  as.data.frame(as.matrix(mcmc_obj)) %>%
    mutate(.iter = row_number())
}

df_int <- mcmc_to_df(res_integrated$samples) %>% mutate(method = "integrated")
df_surv <- mcmc_to_df(res_surv$samples) %>% mutate(method = "survival_only")
df_conn <- mcmc_to_df(res_conn$samples) %>% mutate(method = "connectivity_only")

# On va extraire les paramètres d'intérêt : phi[k], p, mu[k,j]
# normaliser les noms de colonnes pour retrouver phi[1], mu[1,1] etc.

# Fonction pratique pour pivot longer des parametres scalar/matrix nommés par nimble
gather_params <- function(df, param_pattern, param_label) {
  cols <- grep(param_pattern, names(df), value = TRUE)
  if (length(cols) == 0) return(NULL)
  df %>%
    select(all_of(cols), .iter, method) %>%
    pivot_longer(cols = -c(.iter, method), names_to = "param", values_to = "value") %>%
    mutate(type = param_label)
}

# Extractions
phis <- bind_rows(
  gather_params(df_int, "^phi\\[", "phi"),
  gather_params(df_surv, "^phi\\[", "phi"),
  # conn-only has no phi
  NULL
)

ps <- bind_rows(
  gather_params(df_int, "^p$", "p"),
  gather_params(df_surv, "^p$", "p")
)

mus <- bind_rows(
  gather_params(df_int, "^mu\\[", "mu"),
  gather_params(df_conn, "^mu\\[", "mu")
)

sigmas <- bind_rows(
  gather_params(df_int, "^sigma\\[", "sigma"),
  gather_params(df_conn, "^sigma\\[", "sigma")
)

pis <- bind_rows(
  gather_params(df_int, "^pi\\[", "pi"),
  gather_params(df_surv, "^pi\\[", "pi"),
  gather_params(df_conn, "^pi\\[", "pi")
)

# parse indices : mu[1,2] -> cluster 1, dim 2
parse_indices <- function(param_string) {
  # extract numbers between brackets
  nums <- regmatches(param_string, gregexpr("\\d+", param_string))[[1]] %>% as.integer()
  return(nums)
}

# ajouter colonnes utiles
process_mu <- function(df_mu) {
  df_mu %>%
    mutate(nums = map(param, parse_indices)) %>%
    mutate(cluster = map_int(nums, 1),
           dim = map_int(nums, 2)) %>%
    select(-nums)
}
process_phi <- function(df_phi) {
  df_phi %>%
    mutate(cluster = as.integer(gsub("phi\\[|\\]", "", param)))
}
process_sigma <- function(df_sigma) {
  df_sigma %>%
    mutate(nums = map(param, parse_indices)) %>%
    mutate(cluster = map_int(nums, 1),
           dim = map_int(nums, 2)) %>%
    select(-nums)
}
process_pi <- function(df_pi) {
  df_pi %>% mutate(cluster = as.integer(gsub("pi\\[|\\]", "", param)))
}

mus2 <- process_mu(mus)
phis2 <- process_phi(phis)
sigmas2 <- process_sigma(sigmas)
pis2 <- process_pi(pis)

# --------------------------
# 7) RÉSUMÉS POSTÉRIEURS (moyenne, IC 95%) et métriques (bias, RMSE, largeur IC)
# --------------------------

# Function to compute summary per group (by method + param specifics)
summarize_param <- function(df, true_value, group_cols) {
  df %>%
    group_by(across(all_of(group_cols)), method) %>%
    summarise(mean = mean(value),
              sd = sd(value),
              lower = quantile(value, 0.025),
              upper = quantile(value, 0.975),
              width = upper - lower,
              .groups = "drop") %>%
    mutate(true = true_value,
           bias = mean - true,
           rmse = sqrt((mean - true)^2 + sd^2))
}

# PHI summaries: true phi per cluster
phi_summaries <- phis2 %>%
  group_by(method, cluster) %>%
  summarise(mean = mean(value), sd = sd(value),
            lower = quantile(value, .025), upper = quantile(value, .975),
            width = upper - lower, .groups = "drop") %>%
  mutate(true = phi_true[cluster],
         bias = mean - true,
         rmse = sqrt((mean - true)^2 + sd^2))

# MU summaries: true mu per cluster & dim
mu_summaries <- mus2 %>%
  group_by(method, cluster, dim) %>%
  summarise(mean = mean(value), sd = sd(value),
            lower = quantile(value, .025), upper = quantile(value, .975),
            width = upper - lower, .groups = "drop") %>%
  mutate(true = map2_dbl(cluster, dim, ~ mu_true[[.x]][.y]),
         bias = mean - true,
         rmse = sqrt((mean - true)^2 + sd^2))

# SIGMA summaries: true sigma diag per cluster/dim
sigma_true_vec <- unlist(lapply(Sigma_true, diag))
# but Sigma_true is a list of matrices; for matching we'll set true sigma as diag of each Sigma_true
sigma_summaries <- sigmas2 %>%
  group_by(method, cluster, dim) %>%
  summarise(mean = mean(value), sd = sd(value),
            lower = quantile(value, .025), upper = quantile(value, .975),
            width = upper - lower, .groups = "drop") %>%
  mutate(true = map2_dbl(cluster, dim, ~ Sigma_true[[.x]][.y, .y]),
         bias = mean - true,
         rmse = sqrt((mean - true)^2 + sd^2))

# --------------------------
# 8) GRAPHIQUES
# --------------------------

# Helper pour palette
method_colors <- c("integrated" = "#1b9e77", "survival_only" = "#d95f02", "connectivity_only" = "#7570b3")

# 8.a) Densités postérieures phi (overlay), une facette par cluster
plot_phi_density <- function(phis2) {
  ggplot(phis2, aes(x = value, color = method, fill = method)) +
    geom_density(alpha = 0.2) +
    facet_wrap(~ cluster, scales = "free") +
    geom_vline(data = data.frame(cluster = 1:K, true = phi_true), aes(xintercept = true), linetype = "dashed") +
    scale_color_manual(values = method_colors) +
    scale_fill_manual(values = method_colors) +
    labs(title = "Densités postérieures de phi par méthode et cluster",
         x = expression(phi), y = "Densité") +
    theme_minimal()
}

p_phi_density <- plot_phi_density(phis2)

# 8.b) Densités postérieures pour mu (on montre une dimension par facette)
# Nous allons montrer mu dimension 1..d pour chaque cluster (facettage  by cluster x dim)
plot_mu_density <- function(mus2) {
  ggplot(mus2, aes(x = value, color = method, fill = method)) +
    geom_density(alpha = 0.2) +
    facet_grid(cluster ~ dim, scales = "free") +
    geom_vline(data = mu_summaries, aes(xintercept = true), linetype = "dashed") +
    scale_color_manual(values = method_colors) +
    scale_fill_manual(values = method_colors) +
    labs(title = "Densités postérieures de mu (toutes dimensions & clusters)",
         x = "mu", y = "Densité") +
    theme_minimal()
}

p_mu_density <- plot_mu_density(mus2)

# 8.c) Forest plot (IC 95%) pour phi par méthode
plot_phi_forest <- function(phi_summaries) {
  phi_summaries %>%
    ggplot(aes(x = factor(cluster), y = mean, ymin = lower, ymax = upper, color = method)) +
    geom_pointrange(position = position_dodge(width = 0.6)) +
    geom_hline(aes(yintercept = true), data = data.frame(cluster = 1:K, true = phi_true), linetype = "dashed") +
    scale_color_manual(values = method_colors) +
    labs(title = "Estimation de phi : moyennes post. et IC 95% par méthode",
         x = "Cluster", y = expression(phi)) +
    theme_minimal()
}

p_phi_forest <- plot_phi_forest(phi_summaries)

# 8.d) Boxplot des widths (précision) pour mu et phi
plot_widths <- function(mu_summaries, phi_summaries) {
  mu_w <- mu_summaries %>% mutate(param = paste0("mu_c", cluster, "_d", dim)) %>%
    select(method, param, width)
  phi_w <- phi_summaries %>% mutate(param = paste0("phi_c", cluster)) %>%
    select(method, param, width)
  widths <- bind_rows(mu_w, phi_w)
  ggplot(widths, aes(x = method, y = width, fill = method)) +
    geom_boxplot() +
    scale_fill_manual(values = method_colors) +
    labs(title = "Largeur des IC95% (indicateur de précision) — mu & phi",
         x = "Méthode", y = "Largeur IC95%") +
    theme_minimal()
}

p_widths <- plot_widths(mu_summaries, phi_summaries)

# 8.e) Table résumée pour phi : bias & RMSE
phi_metrics_table <- phi_summaries %>%
  select(method, cluster, true, mean, bias, rmse, width) %>%
  arrange(cluster, method)

# 8.f) Traceplots pour phi (premiers 3 clusters)
plot_trace_phi <- function(res_list) {
  # combine sample matrices with method label
  df_i <- as.matrix(res_integrated$samples)[, grep("^phi\\[", colnames(as.matrix(res_integrated$samples))), drop = FALSE]
  df_s <- as.matrix(res_surv$samples)[, grep("^phi\\[", colnames(as.matrix(res_surv$samples))), drop = FALSE]
  # conn has none
  df_i <- as.data.frame(df_i); df_i$iter <- 1:nrow(df_i); df_i$method <- "integrated"
  df_s <- as.data.frame(df_s); df_s$iter <- 1:nrow(df_s); df_s$method <- "survival_only"
  long_i <- df_i %>% pivot_longer(cols = starts_with("phi"), names_to = "param", values_to = "value")
  long_s <- df_s %>% pivot_longer(cols = starts_with("phi"), names_to = "param", values_to = "value")
  long <- bind_rows(long_i, long_s)
  long %>% mutate(cluster = as.integer(gsub("phi\\[|\\]", "", param))) %>%
    ggplot(aes(x = iter, y = value, color = method)) +
    geom_line(alpha = 0.8) +
    facet_wrap(~ cluster, scales = "free_y") +
    scale_color_manual(values = method_colors) +
    labs(title = "Traceplots de phi (integrated vs survival_only)", x = "Itération", y = "phi") +
    theme_minimal()
}

p_phi_trace <- plot_trace_phi(list(res_integrated, res_surv))

# --------------------------
# 9) AFFICHAGE / EXPORT
# --------------------------
# Afficher les plots dans la session
print(p_phi_density)
print(p_mu_density)
print(p_phi_forest)
print(p_widths)
print(p_phi_trace)

# Sauvegarder en fichiers
ggsave("phi_density.png", p_phi_density, width = 10, height = 4)
ggsave("mu_density.png", p_mu_density, width = 10, height = 8)
ggsave("phi_forest.png", p_phi_forest, width = 8, height = 4)
ggsave("widths_ic95.png", p_widths, width = 6, height = 4)
ggsave("phi_trace.png", p_phi_trace, width = 10, height = 4)

# Afficher table résumé phi
print("Résumé phi par méthode (mean, true, bias, rmse, width):")
print(phi_metrics_table)

# Exporter résumés mu en csv
write.csv(mu_summaries, "mu_summaries.csv", row.names = FALSE)
write.csv(phi_summaries, "phi_summaries.csv", row.names = FALSE)
write.csv(sigma_summaries, "sigma_summaries.csv", row.names = FALSE)

# --------------------------
# 10) INTERPRÉTATION RAPIDE (print)
# --------------------------
cat("\nInterprétation rapide (à vérifier visuellement dans les graphiques) :\n")
cat("- Compare les densités & les lignes pointillées (valeurs vraies) : tu devrais voir peu ou pas de biais pour les méthodes.\n")
cat("- Compare la largeur des IC95% : la méthode 'integrated' devrait généralement montrer des IC plus étroits (meilleure précision) pour phi et mu comparé aux méthodes séparées.\n")
cat("- Regarde aussi la variance postérieure (sd) et RMSE dans les fichiers csv pour un bilan numérique.\n\n")

cat("Terminé. Figures sauvegardées : phi_density.png, mu_density.png, phi_forest.png, widths_ic95.png, phi_trace.png\n")
