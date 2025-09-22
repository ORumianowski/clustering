# --------------------------------------------------
# 0) Packages
# --------------------------------------------------
library(nimble)
library(coda)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(cowplot)

set.seed(123)

# --------------------------------------------------
# 1) Données
# --------------------------------------------------
nb_clusters <- 3

load("dataset/soft_connectivity.data.RData")

X <- cbind(
  soft_connectivity.data$longitude_breed,
  soft_connectivity.data$latitude_breed,
  soft_connectivity.data$longitude_wint,
  soft_connectivity.data$latitude_wint
) |> as.matrix()

colnames(X) <- c("lon_breed","lat_breed","lon_wint","lat_wint")

N <- nrow(X)
d <- ncol(X)
K <- nb_clusters

cat("Nombre de NA dans X :", sum(is.na(X)), "\n")

# --------------------------------------------------
# 2) Modèle NIMBLE
# --------------------------------------------------
code_fullcov <- nimbleCode({
  pi[1:K] ~ ddirch(alpha[1:K])
  for (k in 1:K) {
    mu[k,1:d] ~ dmnorm(mu0[1:d], prec = Prec0[1:d,1:d])
    Omega[k,1:d,1:d] ~ dwish(R[1:d,1:d], df)
    Sigma[k,1:d,1:d] <- inverse(Omega[k,1:d,1:d])
  }
  for (i in 1:N) {
    cl[i] ~ dcat(pi[1:K])
    X[i,1:d] ~ dmnorm(mu[cl[i],1:d], prec = Omega[cl[i],1:d,1:d])
  }
})

# --------------------------------------------------
# 3) Constantes & inits
# --------------------------------------------------
xbar <- colMeans(X, na.rm = TRUE)
xvar <- apply(X, 2, var, na.rm = TRUE)

mu0   <- as.numeric(xbar)
Prec0 <- diag(1 / (10^2), d)

R     <- diag(1 / pmax(xvar, 1e-6), d)
df    <- d + 2

constants <- list(K = K, N = N, d = d,
                  alpha = rep(1, K),
                  mu0 = mu0, Prec0 = Prec0,
                  R = R, df = df)

set.seed(123)
init_mu <- matrix(
  rep(xbar, each = K) + matrix(rnorm(K * d, 0, sd = 0.5), nrow = K),
  nrow = K, byrow = FALSE
)

Omega_init <- array(0, dim = c(K, d, d))
for (k in 1:K) {
  Omega_init[k,,] <- diag(1 / apply(X, 2, var, na.rm = TRUE))
}

inits <- list(
  cl   = sample(1:K, N, replace = TRUE),
  mu   = init_mu,
  Omega = Omega_init,
  pi   = rep(1/K, K),
  X = X
)

# --------------------------------------------------
# 4) Helper run_nimble
# --------------------------------------------------
run_nimble <- function(code, constants, inits, monitors,
                       niter = 6000, nburnin = 2000, thin = 5, progressBar = TRUE) {
  model  <- nimbleModel(code, constants = constants, inits = inits, calculate = TRUE)
  cmodel <- compileNimble(model)
  conf   <- configureMCMC(model, monitors = monitors)
  mcmc   <- buildMCMC(conf)
  cmcmc  <- compileNimble(mcmc, project = model)
  samp   <- runMCMC(cmcmc, niter = niter, nburnin = nburnin,
                    thin = thin, progressBar = progressBar)
  list(model = model, cmodel = cmodel, samples = as.mcmc(samp))
}

# --------------------------------------------------
# 5) Lancer MCMC
# --------------------------------------------------
monitors <- c("pi", "mu", "Omega")
res <- run_nimble(code_fullcov, constants, inits, monitors)

samples_df <- as.data.frame(as.matrix(res$samples))

# --------------------------------------------------
# 6) Extraction mu et Sigma
# --------------------------------------------------
mu_cols <- grep("^mu\\[", colnames(samples_df), value = TRUE)
k_idx <- as.integer(sub("^mu\\[(\\d+),\\s*(\\d+)\\]$", "\\1", mu_cols))
d_idx <- as.integer(sub("^mu\\[(\\d+),\\s*(\\d+)\\]$", "\\2", mu_cols))

mu_post <- matrix(NA_real_, nrow = K, ncol = d)
for (c in seq_along(mu_cols)) {
  mu_post[k_idx[c], d_idx[c]] <- mean(samples_df[[mu_cols[c]]])
}
colnames(mu_post) <- colnames(X)

Omega_means <- vector("list", K)
for (k in 1:K) {
  Omega_k <- matrix(NA_real_, d, d)
  for (i in 1:d) {
    for (j in 1:d) {
      colname <- paste0("Omega[", k, ", ", i, ", ", j, "]")
      if (colname %in% colnames(samples_df)) {
        Omega_k[i,j] <- mean(samples_df[[colname]])
      }
    }
  }
  Omega_k <- (Omega_k + t(Omega_k)) / 2
  Omega_means[[k]] <- Omega_k
}

make_spd <- function(S) {
  S <- (S + t(S)) / 2
  ev <- eigen(S)
  vals <- pmax(ev$values, 1e-8)
  S2 <- ev$vectors %*% diag(vals) %*% t(ev$vectors)
  (S2 + t(S2)) / 2
}

Sigma_post <- lapply(Omega_means, function(Om) {
  Om <- make_spd(Om)
  Si <- solve(Om)
  make_spd(Si)
})

# --------------------------------------------------
# 7) Fonctions pour ellipses
# --------------------------------------------------
ellipse_points <- function(mu, Sigma, r, n = 200) {
  ang <- seq(0, 2*pi, length.out = n)
  circle <- cbind(cos(ang), sin(ang))
  R <- chol(Sigma)
  pts <- circle %*% t(R) * r
  sweep(pts, 2, mu, FUN = "+")
}

make_aplats <- function(mu_list, Sigma_list,
                        probs = c(0.60, 0.75, 0.85, 0.90, 0.95),
                        alphas = seq(0.25, 0.15, length.out = 5),
                        fill_cols) {
  r_levels <- sqrt(qchisq(probs, df = 2))
  out <- list(); idx <- 1L
  for (k in seq_along(mu_list)) {
    for (i in seq_along(r_levels)) {
      pts <- ellipse_points(mu_list[[k]], Sigma_list[[k]], r_levels[i])
      out[[idx]] <- data.frame(
        x = pts[,1], y = pts[,2],
        comp = factor(k),
        level = i,
        alpha = alphas[i],
        fill = fill_cols[k]
      )
      idx <- idx + 1L
    }
  }
  df <- bind_rows(out)
  arrange(df, desc(level))
}

# --------------------------------------------------
# 8) Deux cartes Europe (été/hiver)
# --------------------------------------------------
cluster_cols <- c("red3", "forestgreen", "royalblue")

# Été
mu_summer    <- mu_post[, 1:2, drop = FALSE]
Sigma_summer <- lapply(Sigma_post, function(S) S[1:2, 1:2])
polys_summer <- make_aplats(split(mu_summer, row(mu_summer)),
                            Sigma_summer, fill_cols = cluster_cols)
df_cent_summer <- data.frame(x = mu_summer[,1], y = mu_summer[,2], comp = factor(1:K))

# Hiver
mu_winter    <- mu_post[, 3:4, drop = FALSE]
Sigma_winter <- lapply(Sigma_post, function(S) S[3:4, 3:4])
polys_winter <- make_aplats(split(mu_winter, row(mu_winter)),
                            Sigma_winter, fill_cols = cluster_cols)
df_cent_winter <- data.frame(x = mu_winter[,1], y = mu_winter[,2], comp = factor(1:K))

world <- ne_countries(scale = "medium", returnclass = "sf")
base_map <- geom_sf(data = world, fill = "grey92", color = "grey70", size = 0.3)

p_summer <- ggplot() +
  base_map +
  geom_polygon(data = polys_summer,
               aes(x = x, y = y, group = interaction(comp, level),
                   alpha = alpha, fill = fill), color = NA) +
  scale_alpha_identity() +
  scale_fill_identity() +
  geom_point(data = df_cent_summer, aes(x = x, y = y, color = comp),
             shape = 18, size = 3) +
  scale_color_manual(values = cluster_cols, name = "component",
                     labels = paste("Cluster", 1:K)) +
  coord_sf(xlim = c(-15, 35), ylim = c(35, 70), expand = FALSE) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right") +
  labs(title = "Barycentres et ellipses - Été", x = "Longitude", y = "Latitude")

p_winter <- ggplot() +
  base_map +
  geom_polygon(data = polys_winter,
               aes(x = x, y = y, group = interaction(comp, level),
                   alpha = alpha, fill = fill), color = NA) +
  scale_alpha_identity() +
  scale_fill_identity() +
  geom_point(data = df_cent_winter, aes(x = x, y = y, color = comp),
             shape = 18, size = 3) +
  scale_color_manual(values = cluster_cols, name = "component",
                     labels = paste("Cluster", 1:K)) +
  coord_sf(xlim = c(-15, 35), ylim = c(35, 70), expand = FALSE) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "right") +
  labs(title = "Barycentres et ellipses - Hiver", x = "Longitude", y = "Latitude")

plot_grid(p_summer, p_winter, ncol = 2)
