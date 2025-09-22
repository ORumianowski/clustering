# ==========================================================
# 0) Librairies
# ==========================================================
library(nimble)
library(MASS)
library(coda)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(cowplot)
library(gtools)
library(RColorBrewer)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)

set.seed(123)

# ==========================================================
# 1) Jeu de données
# ==========================================================
nb_clusters <- 4   # <-- tu peux changer ici

load("dataset/soft_connectivity.data.RData")

gps <- read.csv("dataset/curlew_metadata_22-02-2024_coordinates.csv", sep = ",", dec = ".", header = TRUE)



X_cmr <- cbind(
  soft_connectivity.data$longitude_breed,
  soft_connectivity.data$latitude_breed,
  soft_connectivity.data$longitude_wint,
  soft_connectivity.data$latitude_wint
) |> as.data.frame()|> 
  na.omit()|> 
  mutate(source = "cmr")


X_gps <- cbind(
  gps$Breeding.site.longitude,
  gps$Breeding.site.latitude,
  gps$Wintering.site.1.longitude,
  gps$Wintering.site.1.latitude
) |> as.data.frame() |> 
  na.omit() |> 
  mutate(source = "gps")

X_full <- rbind(X_cmr,
        X_gps)
X = X_full[,1:4]

# Tableau uniquement pour stats descriptives
X_tbl <- as.data.frame(X)

# ==========================================================
# 2) Modèle NIMBLE
# ==========================================================
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
    for(j in 1:d){
      X[i,j] ~ dnorm(mu[cl[i],j], prec[cl[i],j])
    }
  }
})

# ==========================================================
# 3) Constantes, données, inits
# ==========================================================
N <- nrow(X)
d <- ncol(X)
K <- nb_clusters

# Moyennes des colonnes (pour initialiser mu)
d_means <- colMeans(X_tbl, na.rm = TRUE)

# On génère nb_clusters vecteurs mu autour des moyennes observées
init_mu <- replicate(K, d_means + rnorm(d, 0, 1), simplify = TRUE)
init_mu <- t(init_mu)   # matrice K x d

constants <- list(
  K = K,
  N = N,
  d = d,
  alpha = rep(1, K)
)

data_conn_only <- list(X = X)

inits_conn <- list(
  cl = sample(1:K, N, replace = TRUE),
  mu = init_mu,
  sigma = matrix(1, nrow = K, ncol = d),
  pi = rep(1/K, K)
)

# ==========================================================
# 4) Helper pour lancer Nimble
# ==========================================================
run_nimble <- function(code, constants, data, inits, monitors,
                       niter = 6000, nburnin = 2000, thin = 4,
                       progressBar = TRUE) {
  model  <- nimbleModel(code, constants = constants,
                        data = data, inits = inits, calculate = TRUE)
  cmodel <- compileNimble(model)
  conf   <- configureMCMC(model, monitors = monitors)
  mcmc   <- buildMCMC(conf)
  cmcmc  <- compileNimble(mcmc, project = model)
  samp   <- runMCMC(cmcmc, niter = niter, nburnin = nburnin,
                    thin = thin, progressBar = progressBar)
  list(model = model, cmodel = cmodel, samples = as.mcmc(samp))
}

# ==========================================================
# 5) Lancer le modèle
# ==========================================================
cat("Running connectivity-only model...\n")

res_conn <- run_nimble(
  code_connectivity_only,
  constants,
  data_conn_only,
  inits_conn,
  monitors = c("pi", "mu", "sigma"),
  niter = 6000, nburnin = 2000, thin = 4
)

# --------------------------
# 6) Convertir samples en data.frame
# --------------------------
samples_conn_df <- as.data.frame(as.matrix(res_conn$samples))
samples_conn_df[1:3,]

# --------------------------------------------------
# 1) Extraire mu_post et Sigma_post (diagonaux)
# --------------------------------------------------
mu_cols <- grep("^mu\\[", colnames(samples_conn_df), value = TRUE)
k_idx <- as.integer(sub("^mu\\[(\\d+),\\s*(\\d+)\\]$", "\\1", mu_cols))
d_idx <- as.integer(sub("^mu\\[(\\d+),\\s*(\\d+)\\]$", "\\2", mu_cols))

K <- max(k_idx); D <- max(d_idx)

mu_post <- matrix(NA_real_, nrow = K, ncol = D)
for (c in seq_along(mu_cols)) {
  mu_post[k_idx[c], d_idx[c]] <- mean(samples_conn_df[[mu_cols[c]]])
}
colnames(mu_post) <- c("lon_breed","lat_breed","lon_wint","lat_wint")

# matrices diagonales
Sigma_post <- vector("list", K)
for (k in 1:K) {
  sigmas <- numeric(D)
  for (d in 1:D) {
    colname <- paste0("sigma[", k, ", ", d, "]")
    sigmas[d] <- mean(samples_conn_df[[colname]])
  }
  Sigma_post[[k]] <- diag(sigmas, D, D)
}

# --------------------------------------------------
# 2) Fonctions ellipses/aplats
# --------------------------------------------------
ellipse_points <- function(mu, Sigma, r, n = 200) {
  ang <- seq(0, 2*pi, length.out = n)
  circle <- cbind(cos(ang), sin(ang))
  R <- chol(Sigma)
  pts <- circle %*% t(R) * r
  sweep(pts, 2, mu, FUN = "+")
}

make_aplats <- function(mu_list, Sigma_list,
                        probs = c(0.50, 0.75, 0.90, 0.95, 0.99),
                        alphas = seq(0.15, 0.15, length.out = length(probs)),
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
  dplyr::bind_rows(out) |> dplyr::arrange(desc(level))
}

# Données brutes pour affichage
df_points_summer <- X_tbl |> 
  dplyr::select(lon_breed = V1, lat_breed = V2)

df_points_winter <- X_tbl |> 
  dplyr::select(lon_wint = V3, lat_wint = V4)

# Palette de couleurs dynamique
cluster_cols <- brewer.pal(max(3, K), "Set2")[1:K]

# Aplats
mu_summer <- mu_post[,1:2, drop=FALSE]
mu_winter <- mu_post[,3:4, drop=FALSE]

Sigma_summer <- lapply(Sigma_post, function(M) M[1:2,1:2])
Sigma_winter <- lapply(Sigma_post, function(M) M[3:4,3:4])

polys_summer <- make_aplats(split(mu_summer, row(mu_summer)),
                            Sigma_summer, fill_cols = cluster_cols)

polys_winter <- make_aplats(split(mu_winter, row(mu_winter)),
                            Sigma_winter, fill_cols = cluster_cols)

# Centroides
df_cent_summer <- data.frame(
  x = mu_post[,1], y = mu_post[,2],
  comp = factor(1:K)
)
df_cent_winter <- data.frame(
  x = mu_post[,3], y = mu_post[,4],
  comp = factor(1:K)
)

# ==========================================================
# Base map
# ==========================================================
world <- ne_countries(scale = "medium", returnclass = "sf")
base_map <- geom_sf(data = world, fill = "antiquewhite")

# ==========================================================
# Plot final
# ==========================================================
# Été (Breeding)
p_summer <- ggplot() +
  base_map +
  geom_point(data = df_points_summer,
             aes(x = lon_breed, y = lat_breed),
             color = "black", alpha = 0.5, size = 1.5) +
  geom_polygon(data = polys_summer,
               aes(x = x, y = y, group = interaction(comp, level),
                   alpha = alpha, fill = fill), color = NA) +
  scale_alpha_identity() +
  scale_fill_identity() +
  geom_point(data = df_cent_summer, aes(x = x, y = y, color = comp),
             shape = 18, size = 3) +
  scale_color_manual(values = cluster_cols,
                     labels = paste("Cluster", 1:K), name = NULL) +
  coord_sf(xlim = c(-15, 35), ylim = c(35, 70), expand = FALSE) +
  theme_minimal(base_size = 13) +
  labs(title = "Breeding",
       x = "Longitude", y = "Latitude")

# Hiver (Wintering)
p_winter <- ggplot() +
  base_map +
  geom_point(data = df_points_winter,
             aes(x = lon_wint, y = lat_wint),
             color = "black", alpha = 0.5, size = 1.5) +
  geom_polygon(data = polys_winter,
               aes(x = x, y = y, group = interaction(comp, level),
                   alpha = alpha, fill = fill), color = NA) +
  scale_alpha_identity() +
  scale_fill_identity() +
  geom_point(data = df_cent_winter, aes(x = x, y = y, color = comp),
             shape = 18, size = 3) +
  scale_color_manual(values = cluster_cols,
                     labels = paste("Cluster", 1:K), name = NULL) +
  coord_sf(xlim = c(-15, 35), ylim = c(35, 70), expand = FALSE) +
  theme_minimal(base_size = 13) +
  labs(title = "Wintering",
       x = "Longitude", y = "Latitude")

# ==========================================================
# Affichage côte à côte
# ==========================================================
final_plot <- plot_grid(p_summer, p_winter, ncol = 2)

# Affichage à l'écran
print(final_plot)

# ==========================================================
# Sauvegarde en haute qualité
# ==========================================================
ggsave(
  filename = paste0("clusters_", K, "_map.png"),
  plot = final_plot,
  width = 14, height = 7, dpi = 300
)
