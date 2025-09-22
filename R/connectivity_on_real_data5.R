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
nb_clusters <- 3   

load("dataset/soft_connectivity.data.RData")

gps <- read.csv("dataset/curlew_metadata_22-02-2024_coordinates.csv",
                sep = ",", dec = ".", header = TRUE)

X_cmr <- cbind(
  soft_connectivity.data$longitude_breed,
  soft_connectivity.data$latitude_breed,
  soft_connectivity.data$longitude_wint,
  soft_connectivity.data$latitude_wint
) |> as.data.frame() |>
  na.omit() |>
  mutate(source = "cmr")

X_gps <- cbind(
  gps$Breeding.site.longitude,
  gps$Breeding.site.latitude,
  gps$Wintering.site.1.longitude,
  gps$Wintering.site.1.latitude
) |> as.data.frame() |>
  na.omit() |>
  mutate(source = "gps")

X_full <- rbind(X_cmr, X_gps)
X <- X_full[,1:4]

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

# # --------------------------
# 6) Convertir samples en data.frame
# --------------------------
samples_conn_df <- as.data.frame(as.matrix(res_conn$samples))

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

# Données brutes pour affichage
df_points_summer <- X_full |> 
  dplyr::select(lon_breed = V1, lat_breed = V2, source)

df_points_winter <- X_full |> 
  dplyr::select(lon_wint = V3, lat_wint = V4, source)

# Palette de couleurs clusters
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
# ==========================================================
# 0) Librairies
# ==========================================================
library(dplyr)
library(ggplot2)
library(viridisLite)
library(patchwork)
library(rnaturalearth)
library(rnaturalearthdata)

# ==========================================================
# 1) Palette et constantes
# ==========================================================
cluster_levels  <- factor(seq_len(K), levels = as.character(seq_len(K)))
cluster_palette <- setNames(viridisLite::viridis(K, option = "D"),
                            levels(cluster_levels))

coverage_levels <- c("50%", "75%", "90%", "95%", "99%")
coverage_alphas <- c("50%" = 0.25,
                     "75%" = 0.15,
                     "90%" = 0.10,
                     "95%" = 0.08,
                     "99%" = 0.05)
# ==========================================================
# 2) Fonctions ellipses
# ==========================================================
ellipse_points <- function(mu, Sigma, r, n = 200) {
  ang <- seq(0, 2*pi, length.out = n)
  circle <- cbind(cos(ang), sin(ang))
  R <- chol(Sigma)
  pts <- circle %*% t(R) * r
  sweep(pts, 2, mu, FUN = "+")
}

make_aplats <- function(mu_list, Sigma_list,
                        probs = c(0.50, 0.75, 0.90, 0.95, 0.99)) {
  r_levels  <- sqrt(qchisq(probs, df = 2))
  cov_labels <- sprintf("%d%%", round(100*probs))
  out <- list(); idx <- 1L
  for (k in seq_along(mu_list)) {
    for (i in seq_along(r_levels)) {
      pts <- ellipse_points(mu_list[[k]], Sigma_list[[k]], r_levels[i])
      out[[idx]] <- data.frame(
        x        = pts[,1],
        y        = pts[,2],
        cluster  = factor(k, levels = levels(cluster_levels)),
        coverage = factor(cov_labels[i], levels = coverage_levels)
      )
      idx <- idx + 1L
    }
  }
  bind_rows(out)
}

# ==========================================================
# 3) Fonction de construction de carte
# ==========================================================

# ==========================================================
# Thème publication avec Helvetica
# ==========================================================
theme_pub <- theme_minimal(base_size = 12) +
  theme(
    text             = element_text(family = "Helvetica"),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.title       = element_text(face = "bold"),
    axis.text        = element_text(color = "grey10"),
    plot.title       = element_text(face = "bold", hjust = 0.5),
    legend.position  = "bottom",
    legend.box       = "vertical",
    legend.title     = element_text(face = "bold")
  )

build_map <- function(points_df, polys_df, cents_df, title) {
  ggplot() +
    base_map +
    
    # --- Points CMR (cercles blancs bordés de noir) ---
    geom_point(data = points_df |> filter(source == "cmr"),
               aes(x = lon, y = lat, shape = "CMR"),
               size = 2, color = "black", fill = "white", stroke = 0.6) +
    
    # --- Points GPS (losanges blancs bordés de noir) ---
    geom_point(data = points_df |> filter(source == "gps"),
               aes(x = lon, y = lat, shape = "GPS"),
               size = 2, color = "black", fill = "white", stroke = 0.6) +
    
    # --- Ellipses translucides ---
    geom_polygon(
      data = polys_df,
      aes(x = x, y = y,
          group = interaction(cluster, coverage),
          fill  = cluster,
          alpha = coverage),
      color = NA
    ) +
    scale_fill_manual(values = cluster_palette, name = "Clusters") +
    scale_alpha_manual(values = coverage_alphas,
                       name = "Probability contours") +
    
    # --- Centroïdes ---
    geom_point(data = cents_df, aes(x = x, y = y, color = cluster),
               shape = 3, size = 2.5, stroke = 1, show.legend = FALSE) +
    scale_color_manual(values = cluster_palette, guide = "none") +
    
    # --- Légende des sources ---
    scale_shape_manual(values = c("CMR" = 21, "GPS" = 23),
                       name = "Data source") +
    
    guides(
      shape = guide_legend(order = 1,
                           override.aes = list(fill = c("white","white"),
                                               color = "black", size = 2.8, stroke = 0.8)),
      fill  = guide_legend(order = 2),
      alpha = guide_legend(order = 3)
    ) +
    
    coord_sf(xlim = c(-15, 35), ylim = c(35, 70), expand = FALSE) +
    labs(title = title, x = "Longitude", y = "Latitude") +
    theme_pub
}

# ==========================================================
# 4) Cartes et assemblage
# ==========================================================
p_summer <- build_map(df_points_summer, polys_summer, df_cent_summer, "Breeding")
p_winter <- build_map(df_points_winter, polys_winter, df_cent_winter, "Wintering")

final_plot <- p_summer + p_winter + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# ==========================================================
# 5) Sauvegarde
# ==========================================================
ggsave("clusters_map_publication.png", final_plot,
       type = "cairo", dpi = 600,
       width = 340, height = 200, units = "mm", bg = "white")



ggsave("clusters_map_publication.png", final_plot,
       type = "cairo", dpi = 600, width = 180, height = 90, units = "mm")
