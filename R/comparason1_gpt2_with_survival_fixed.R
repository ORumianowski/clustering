# compare_integrated_vs_separate_ordered_phi.R
# Script complet : simulation, 3 modèles (integrated, survival-only with ordered phi, connectivity-only),
# realignment of labels, comparisons, plots.
#
# Packages nécessaires :
# install.packages(c("MASS","nimble","coda","ggplot2","dplyr","tidyr","purrr","cowplot","gtools"))
library(nimble)
library(MASS)
library(coda)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(cowplot)
library(gtools)  # pour permutations

set.seed(123)

# --------------------------
# 1) SIMULATION (identique à ton code)
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

phi_true <- c(0.9, 0.7, 0.5)  # attention : pas ordonné nécessairement
p_true <- 0.6

# Génération X
X <- matrix(NA, nrow = N, ncol = d)
for (i in 1:N) {
  k <- z_true[i]
  X[i, ] <- mvrnorm(1, mu_true[[k]], Sigma_true[[k]])
}

# États latents et observations y
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
# 2) CODES NIMBLE
# --------------------------

# 2.a) Intégré (comme avant)
code_integrated <- nimbleCode({
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

# 2.b) Survie seule, pi fixé & phi ordonné
# Paramétrisation pour forcer phi1 < phi2 < phi3 :
# on défini u1,u2,u3 ~ Beta puis on construit phi cumulés croissants
code_survival_ordered <- nimbleCode({
  # pi fixé, on ne l'estime pas ici (on ne la met pas en prior)
  p ~ dbeta(1, 1)
  # increments
  for(k in 1:K) {
    delta[k] ~ dbeta(1,1) # increments dans (0,1)
  }
  # construire phi ordonné croissant entre 0 et 1 :
  phi[1] <- delta[1] * 1.0                        # dans (0,1)
  phi[2] <- phi[1] + delta[2] * (1 - phi[1])      # phi2 > phi1
  phi[3] <- phi[2] + delta[3] * (1 - phi[2])      # phi3 > phi2
  # transitions
  for (k in 1:K) {
    transition[k,1,1] <- phi[k]
    transition[k,1,2] <- 1 - phi[k]
    transition[k,2,1] <- 0
    transition[k,2,2] <- 1
  }
  obs[1,1] <- p
  obs[1,2] <- 1 - p
  obs[2,1] <- 0
  obs[2,2] <- 1
  
  for (i in 1:N) {
    # Ici on fixe z selon pi_true (proportions connues)
    z[i] ~ dcat(pi_fix[1:K])
    state[i,1] ~ dcat(init[1:2])
    for (t in 2:Tmax) {
      state[i,t] ~ dcat(transition[z[i], state[i,t-1], 1:2])
    }
    for (t in 1:Tmax) {
      y[i,t] ~ dcat(obs[state[i,t], 1:2])
    }
  }
})

# 2.c) Connectivité seule
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
data_conn_only <- list(X = X)
# pour modèle survie: fournir pi_fix en tant que data (pi_true fixé)
data_surv_only <- list(y = y, pi_fix = pi_true)

# inits
inits_common <- list(
  z = sample(1:K, N, replace = TRUE),
  mu = matrix(c(0, 0, 0, 0,
                3, 3, 3, 3,
                -3, 3, -3, 3), nrow = K, ncol = d, byrow = TRUE),
  sigma = matrix(rep(1, K * d), nrow = K, ncol = d),
  # pour survie ordonnée : delta inits (dans 0,1)
  delta = rep(0.5, K),
  p = 0.5,
  state = matrix(1, N, Tmax)
)

inits_integrated <- inits_common
inits_conn <- inits_common[c("z", "mu", "sigma")]
inits_conn$pi <- rep(1/K, K)
inits_surv <- list(z = sample(1:K, N, replace = TRUE), delta = rep(0.6, K), p = 0.5, state = matrix(1, N, Tmax))
# (z is used but driven by pi_fix in data)

# --------------------------
# 4) Fonction run nimble
# --------------------------
run_nimble <- function(code, constants, data, inits, monitors,
                       niter = 6000, nburnin = 2000, thin = 4, showProgress = TRUE) {
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
# 5) LANCER LES 3 CHAÎNES
# --------------------------
niter <- 6000
nburnin <- 2000
thin <- 4

cat("Running integrated model...\n")
res_int <- run_nimble(code_integrated, constants, data_integrated, inits_integrated,
                      monitors = c("pi", "mu", "sigma", "phi", "p"),
                      niter = niter, nburnin = nburnin, thin = thin)

cat("Running survival-only (pi fixed, phi ordered)...\n")
res_surv <- run_nimble(code_survival_ordered, constants, data_surv_only, inits_surv,
                       monitors = c("phi", "delta", "p"),
                       niter = niter, nburnin = nburnin, thin = thin)

cat("Running connectivity-only model...\n")
res_conn <- run_nimble(code_connectivity_only, constants, data_conn_only, inits_conn,
                       monitors = c("pi", "mu", "sigma"),
                       niter = niter, nburnin = nburnin, thin = thin)

# --------------------------
# 6) EXTRACTION DES CHAÎNES EN DATA.FRAME
# --------------------------
samples_int_df <- as.data.frame(as.matrix(res_int$samples))
samples_surv_df <- as.data.frame(as.matrix(res_surv$samples))
samples_conn_df <- as.data.frame(as.matrix(res_conn$samples))

# --------------------------
# 7) FONCTIONS D'AIDE : EXTRACTION moyennes/posterior medians et recherche permutation
# --------------------------

# extraire mu estimées (médianes) pour un modèle qui a mu[.,.]
extract_mu_est <- function(samples_df) {
  mu_cols <- grep("^mu\\[", names(samples_df), value = TRUE)
  if(length(mu_cols) == 0) return(NULL)
  # récupérer K x d en prenant la médiane pour chaque mu[k,j]
  # parse names mu[k,j]
  parse_mu_name <- function(name) {
    nums <- regmatches(name, gregexpr("\\d+", name))[[1]] %>% as.integer()
    return(nums)
  }
  # construire matrix
  mu_est <- matrix(NA, nrow = K, ncol = d)
  for(name in mu_cols) {
    nums <- parse_mu_name(name)
    k <- nums[1]; j <- nums[2]
    mu_est[k, j] <- median(samples_df[[name]])
  }
  return(mu_est)
}

# extraire phi estimées (médianes)
extract_phi_est <- function(samples_df) {
  phi_cols <- grep("^phi\\[", names(samples_df), value = TRUE)
  if(length(phi_cols) == 0) return(NULL)
  phi_est <- numeric(K)
  for(name in phi_cols) {
    k <- as.integer(gsub("phi\\[|\\]", "", name))
    phi_est[k] <- median(samples_df[[name]])
  }
  return(phi_est)
}

# trouver meilleure permutation pour rapprocher estims (matrix Kxd) de mu_true_mat (Kxd)
find_best_perm_matrix <- function(est_mat, true_mat) {
  perms <- permutations(n = K, r = K, v = 1:K)
  best_perm <- NULL
  best_err <- Inf
  for(i in seq_len(nrow(perms))) {
    perm <- perms[i, ]
    err <- sum((est_mat[perm, , drop = FALSE] - true_mat)^2)
    if(err < best_err) {
      best_err <- err
      best_perm <- perm
    }
  }
  return(best_perm)
}

# pour vecteurs (phi)
find_best_perm_vector <- function(est_vec, true_vec) {
  perms <- permutations(n = K, r = K, v = 1:K)
  best_perm <- NULL
  best_err <- Inf
  for(i in seq_len(nrow(perms))) {
    perm <- perms[i, ]
    err <- sum((est_vec[perm] - true_vec)^2)
    if(err < best_err) {
      best_err <- err
      best_perm <- perm
    }
  }
  return(best_perm)
}

# --------------------------
# 8) TROUVER LES PERMUTATIONS ET RÉALIGNER LES NOMS DE COLONNE DES SAMPLES
# --------------------------

# construire mu_true matrix
mu_true_mat <- do.call(rbind, mu_true)  # K x d

# integrated: use mu to find permutation
mu_int_est <- extract_mu_est(samples_int_df)
perm_int <- find_best_perm_matrix(mu_int_est, mu_true_mat)
cat("Permutation optimale (integrated):", perm_int, "\n")

# connectivity-only
mu_conn_est <- extract_mu_est(samples_conn_df)
perm_conn <- find_best_perm_matrix(mu_conn_est, mu_true_mat)
cat("Permutation optimale (connectivity-only):", perm_conn, "\n")

# survival-only : no mu; align phi to phi_true (we want phi correspondance to true clusters)
phi_surv_est <- extract_phi_est(samples_surv_df)
perm_surv <- find_best_perm_vector(phi_surv_est, phi_true)
cat("Permutation optimale (survival-only):", perm_surv, "\n")

# fonction pour renommer colonnes dans un dataframe de samples selon permutation perm
rename_samples_by_perm <- function(samples_df, perm) {
  new_names <- names(samples_df)
  # mu[k,j]
  new_names <- sapply(new_names, function(nm) {
    # mu
    if(grepl("^mu\\[", nm)) {
      nums <- regmatches(nm, gregexpr("\\d+", nm))[[1]] %>% as.integer()
      k <- nums[1]; j <- nums[2]
      newk <- perm[k]
      return(paste0("mu[", newk, ",", j, "]"))
    }
    # sigma
    if(grepl("^sigma\\[", nm)) {
      nums <- regmatches(nm, gregexpr("\\d+", nm))[[1]] %>% as.integer()
      k <- nums[1]; j <- nums[2]
      newk <- perm[k]
      return(paste0("sigma[", newk, ",", j, "]"))
    }
    # phi
    if(grepl("^phi\\[", nm)) {
      k <- as.integer(gsub("phi\\[|\\]", "", nm))
      newk <- perm[k]
      return(paste0("phi[", newk, "]"))
    }
    # pi
    if(grepl("^pi\\[", nm)) {
      k <- as.integer(gsub("pi\\[|\\]", "", nm))
      newk <- perm[k]
      return(paste0("pi[", newk, "]"))
    }
    # otherwise keep same
    return(nm)
  }, USE.NAMES = FALSE)
  names(samples_df) <- new_names
  return(samples_df)
}

# appliquer renommage
samples_int_df_aligned <- rename_samples_by_perm(samples_int_df, perm_int)
samples_conn_df_aligned <- rename_samples_by_perm(samples_conn_df, perm_conn)
samples_surv_df_aligned <- rename_samples_by_perm(samples_surv_df, perm_surv)

# --------------------------
# 9) EXTRAIRE PARAMÈTRES D'INTÉRÊT APRÈS ALIGNEMENT (pour graphiques & résumés)
# --------------------------
# Fonction pour récupérer data.frame long de paramètres (param, value, method)
df_from_samples <- function(samples_df, method_label) {
  df <- as.data.frame(samples_df) %>% mutate(.iter = row_number(), method = method_label) %>%
    pivot_longer(cols = -c(.iter, method), names_to = "param", values_to = "value")
  return(df)
}

df_int_long <- df_from_samples(samples_int_df_aligned, "integrated")
df_conn_long <- df_from_samples(samples_conn_df_aligned, "connectivity_only")
df_surv_long <- df_from_samples(samples_surv_df_aligned, "survival_only")

all_long <- bind_rows(df_int_long, df_conn_long, df_surv_long)

# séparer mu, phi, sigma, pi, p
all_long <- all_long %>%
  mutate(type = case_when(
    grepl("^mu\\[", param) ~ "mu",
    grepl("^sigma\\[", param) ~ "sigma",
    grepl("^phi\\[", param) ~ "phi",
    grepl("^pi\\[", param) ~ "pi",
    param == "p" ~ "p",
    TRUE ~ "other"
  ))

# extraire indices pour mu & sigma & phi & pi
extract_indices <- function(param) {
  nums <- regmatches(param, gregexpr("\\d+", param))[[1]] %>% as.integer()
  if(length(nums) == 2) return(list(cluster = nums[1], dim = nums[2]))
  if(length(nums) == 1) return(list(cluster = nums[1], dim = NA_integer_))
  return(list(cluster = NA_integer_, dim = NA_integer_))
}

all_long <- all_long %>%
  rowwise() %>%
  mutate(idx = list(extract_indices(param)),
         cluster = idx$cluster,
         dim = idx$dim) %>%
  ungroup() %>%
  select(-idx)

# --------------------------
# 10) RÉSUMÉS POSTERIEURS (moyenne, IC95%)
# --------------------------
summary_by_param <- all_long %>%
  filter(type %in% c("mu", "phi", "sigma", "pi", "p")) %>%
  group_by(method, type, param, cluster, dim) %>%
  summarise(mean = mean(value),
            sd = sd(value),
            lower = quantile(value, .025),
            upper = quantile(value, .975),
            width = upper - lower,
            .groups = "drop")

# joindre vraies valeurs selon type
# construire table des vraies valeurs pour mu, phi, sigma, pi
true_df_mu <- expand.grid(cluster = 1:K, dim = 1:d) %>%
  mutate(true = map2_dbl(cluster, dim, ~ mu_true[[.x]][.y]),
         type = "mu")
true_df_phi <- data.frame(cluster = 1:K, true = phi_true, type = "phi")
true_df_sigma <- expand.grid(cluster = 1:K, dim = 1:d) %>%
  mutate(true = map2_dbl(cluster, dim, ~ Sigma_true[[.x]][.y, .y]),
         type = "sigma")
true_df_pi <- data.frame(cluster = 1:K, true = pi_true, type = "pi")

true_df <- bind_rows(true_df_mu, true_df_phi, true_df_sigma, true_df_pi)

# merge
summary_with_true <- summary_by_param %>%
  left_join(true_df, by = c("cluster", "type"))

# calcul bias & rmse
summary_with_true <- summary_with_true %>%
  mutate(bias = mean - true,
         rmse = sqrt((mean - true)^2 + sd^2))

# --------------------------
# 11) GRAPHIQUES (densités & forest plots) — avec labels alignés
# --------------------------
method_colors <- c("integrated" = "#1b9e77", "survival_only" = "#d95f02", "connectivity_only" = "#7570b3")

# 11.a) Densité phi par cluster
plot_phi_density <- function(all_long) {
  dat <- all_long %>% filter(type == "phi")
  ggplot(dat, aes(x = value, color = method, fill = method)) +
    geom_density(alpha = 0.25) +
    facet_wrap(~ cluster, scales = "free") +
    geom_vline(data = data.frame(cluster = 1:K, true = phi_true), aes(xintercept = true), linetype = "dashed") +
    scale_color_manual(values = method_colors) + scale_fill_manual(values = method_colors) +
    labs(title = "Densités postérieures de phi (après réalignement)", x = expression(phi), y = "Densité") +
    theme_minimal()
}
p_phi_density <- plot_phi_density(all_long)

# 11.b) Densité mu par cluster x dim
plot_mu_density <- function(all_long) {
  dat <- all_long %>% filter(type == "mu")
  ggplot(dat, aes(x = value, color = method, fill = method)) +
    geom_density(alpha = 0.25) +
    facet_grid(cluster ~ dim, scales = "free") +
    geom_vline(data = true_df_mu, aes(xintercept = true), linetype = "dashed") +
    scale_color_manual(values = method_colors) + scale_fill_manual(values = method_colors) +
    labs(title = "Densités postérieures de mu (après réalignement)", x = "mu", y = "Densité") +
    theme_minimal()
}
p_mu_density <- plot_mu_density(all_long)

# 11.c) Forest plot pour phi
plot_phi_forest <- function(summary_with_true) {
  dat <- summary_with_true %>% filter(type == "phi")
  ggplot(dat, aes(x = factor(cluster), y = mean, ymin = lower, ymax = upper, color = method)) +
    geom_pointrange(position = position_dodge(width = 0.6)) +
    geom_hline(aes(yintercept = true), data = data.frame(cluster = 1:K, true = phi_true), linetype = "dashed") +
    scale_color_manual(values = method_colors) +
    labs(title = "Phi : moyennes post. et IC95% (après réalignement)", x = "Cluster", y = expression(phi)) +
    theme_minimal()
}
p_phi_forest <- plot_phi_forest(summary_with_true)

# 11.d) Largeur IC pour mu & phi (précision)
plot_widths <- function(summary_with_true) {
  dat <- summary_with_true %>% filter(type %in% c("mu", "phi"))
  ggplot(dat, aes(x = method, y = width, fill = method)) +
    geom_boxplot() +
    scale_fill_manual(values = method_colors) +
    labs(title = "Largeur des IC95% (mu & phi) — indicateur de précision", x = "Méthode", y = "Largeur IC95%") +
    theme_minimal()
}
p_widths <- plot_widths(summary_with_true)

# 11.e) Traceplots pour phi (integrated vs survival_only)
plot_phi_trace <- function(samples_int_df_aligned, samples_surv_df_aligned) {
  int_phi_cols <- grep("^phi\\[", names(samples_int_df_aligned), value = TRUE)
  surv_phi_cols <- grep("^phi\\[", names(samples_surv_df_aligned), value = TRUE)
  df_i <- as.data.frame(samples_int_df_aligned)[, int_phi_cols, drop = FALSE] %>% mutate(iter = 1:nrow(.), method = "integrated")
  df_s <- as.data.frame(samples_surv_df_aligned)[, surv_phi_cols, drop = FALSE] %>% mutate(iter = 1:nrow(.), method = "survival_only")
  long_i <- df_i %>% pivot_longer(cols = -c(iter, method), names_to = "param", values_to = "value")
  long_s <- df_s %>% pivot_longer(cols = -c(iter, method), names_to = "param", values_to = "value")
  long <- bind_rows(long_i, long_s)
  long <- long %>% mutate(cluster = as.integer(gsub("phi\\[|\\]", "", param)))
  ggplot(long, aes(x = iter, y = value, color = method)) +
    geom_line(alpha = 0.6) +
    facet_wrap(~ cluster, scales = "free_y") +
    scale_color_manual(values = method_colors) +
    labs(title = "Traceplots de phi (integrated vs survival_only)", x = "Itération", y = "phi") +
    theme_minimal()
}
p_phi_trace <- plot_phi_trace(samples_int_df_aligned, samples_surv_df_aligned)

# --------------------------
# 12) AFFICHAGE ET EXPORT
# --------------------------
print(p_phi_density)
print(p_mu_density)
print(p_phi_forest)
print(p_widths)
print(p_phi_trace)

ggsave("phi_density_aligned.png", p_phi_density, width = 10, height = 4)
ggsave("mu_density_aligned.png", p_mu_density, width = 10, height = 8)
ggsave("phi_forest_aligned.png", p_phi_forest, width = 8, height = 4)
ggsave("widths_ic95_aligned.png", p_widths, width = 6, height = 4)
ggsave("phi_trace_aligned.png", p_phi_trace, width = 10, height = 4)

# Export CSVs
write.csv(summary_with_true, "summary_with_true_aligned.csv", row.names = FALSE)

cat("\nRésumés sauvegardés dans 'summary_with_true_aligned.csv' et figures PNG.\n")
cat("Permutations appliquées:\n")
cat("  integrated perm:", perm_int, "\n")
cat("  connectivity perm:", perm_conn, "\n")
cat("  survival perm:", perm_surv, "\n\n")

# Impression rapide des metrics phi
phi_metrics <- summary_with_true %>% filter(type == "phi") %>% arrange(cluster, method)
print(phi_metrics %>% select(method, cluster, true, mean, bias, rmse, width))

cat("\nFin du script.\n")
