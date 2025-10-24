# ==================================================
# GMM bayésien (NIMBLE) en 4D (été + hiver) avec correction par Ω
# ==================================================

rm(list = ls())
set.seed(123)

suppressPackageStartupMessages({
  library(nimble)
  library(MASS)
  library(Matrix)
  library(terra)
  library(cowplot)
  library(mvtnorm)
  library(ggplot2)
  library(dplyr)
})

# ----------------------- Utilitaires -----------------------
make_cov2D <- function(angle_deg, var_major, var_minor){
  phi <- angle_deg * pi/180
  R <- matrix(c(cos(phi), -sin(phi), sin(phi), cos(phi)), 2, 2, byrow = TRUE)
  S <- R %*% diag(c(var_major, var_minor)) %*% t(R)
  ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  if (!isSymmetric(S) || any(ev <= 1e-10)) S <- as.matrix(Matrix::nearPD(S)$mat)
  S
}

build_sigma <- function(sd, R){
  D4 <- diag(sd, 4, 4)
  S <- D4 %*% R %*% D4
  ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  if (!isSymmetric(S) || any(ev <= 1e-10)) S <- as.matrix(nearPD(S)$mat)
  S
}

gauss2d <- function(x, y, mu, Sigma) {
  invS <- solve(Sigma)
  dx <- x - mu[1]; dy <- y - mu[2]
  Q <- invS[1,1]*dx*dx + (invS[1,2]+invS[2,1])*dx*dy + invS[2,2]*dy*dy
  exp(-0.5 * Q)
}

# ======================= 1) Simulation 4D =======================
# Nombre de clusters
K <- 3
# Dimension
D <- 4
# Taille de la population
N_tot <- 5000
# Répartition de la population entre les clusters
pi_true <- c(0.20, 0.35, 0.45)

# Centres des clusters
# (x_su, y_su, x_wi, y_wi)
mu_true <- list(
  c(-3.5,  0.0,  -1,   2), 
  c( 0.0, -0.6,  -1,   2),  
  c( 3.5,  0.0,  2,  5)  
)
# Décalage global de la coordonnée y_hiver
for (k in 1:K) mu_true[[k]][4] <- mu_true[[k]][4] - 30

# Bloc-diagonal : (1,2) ⟂ (3,4) — covariances été/hiver nulles
sd1 <- c(0.7, 2.8, 1.0, 2.0)
R1 <- matrix(c(
  1.00,  0.00,  0.00,  0.00,
  0.00,  1.00,  0.00,  0.00,
  0.00,  0.00,  1.00,  0.20,
  0.00,  0.00,  0.20,  1.00
), 4, 4, byrow = TRUE)
Sigma1 <- build_sigma(sd1, R1)

sd2 <- c(2.6, 2.4, 1.0, 2.0)
R2 <- matrix(c(
  1.00, -0.70,  0.00,  0.00,
  -0.70,  1.00,  0.00,  0.00,
  0.00,  0.00,  1.00,  0.20,
  0.00,  0.00,  0.20,  1.00
), 4, 4, byrow = TRUE)
Sigma2 <- build_sigma(sd2, R2)

sd3 <- c(0.8, 3.2, 0.6, 2.8)
R3 <- matrix(c(
  1.00,  0.00,  0.00,  0.00,
  0.00,  1.00,  0.00,  0.00,
  0.00,  0.00,  1.00,  0.00,
  0.00,  0.00,  0.00,  1.00
), 4, 4, byrow = TRUE)
Sigma3 <- build_sigma(sd3, R3)

Sigma_true <- list(Sigma1, Sigma2, Sigma3)

# Échantillonnage latent + données 4D
# Attribution du cluster selon pi
z <- sample(1:K, N_tot, replace = TRUE, prob = pi_true)

# Échantillonnage des coordonnées selon la gaussienne du cluster
X4_all <- t(vapply(1:N_tot, function(i) MASS::mvrnorm(1, mu = mu_true[[z[i]]], Sigma = Sigma_true[[z[i]]]), numeric(D)))

colnames(X4_all) <- c("x_su","y_su","x_wi","y_wi")
X_su_all <- X4_all[,1:2, drop = FALSE]
X_wi_all <- X4_all[,3:4, drop = FALSE]

# ======================= 2) Ω été / Ω hiver =======================
# Grilles/rasters séparés pour l'intégration 2D
pad <- 1
xmin_all <- min(c(X_su_all[,1], X_wi_all[,1])) - pad
xmax_all <- max(c(X_su_all[,1], X_wi_all[,1])) + pad
ymin_all <- min(c(X_su_all[,2], X_wi_all[,2])) - pad
ymax_all <- max(c(X_su_all[,2], X_wi_all[,2])) + pad

# Fonctions de construction de l'effort d'échantillonnage Ω (été et hiver)
build_omega_summer <- function(xmin, xmax, ymin, ymax, nx = 60, ny = 60){
  r <- rast(nrows = ny, ncols = nx, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
  xy <- as.data.frame(xyFromCell(r, 1:ncell(r))); colnames(xy) <- c("x","y")
  omega_fun <- function(x, y){
    base <- 0.04
    bumpL <- 0.90 * gauss2d(x, y, c(-3.8, 2.6), diag(c(0.60, 0.60)))
    bumpR <- 0.90 * gauss2d(x, y, c( 3.8, 2.6), diag(c(0.65, 0.65)))
    bumpC <- 1.00 * gauss2d(x, y, c( 0.0,-3.2), diag(c(6.0, 0.9)))
    pmax(0, pmin(base + bumpL + bumpR + bumpC, 0.70))
  }
  r[] <- with(xy, omega_fun(x, y)) * 0.5
  r
}

build_omega_winter <- function(xmin, xmax, ymin, ymax, nx = 60, ny = 60){
  r <- rast(nrows = ny, ncols = nx, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
  xy <- as.data.frame(xyFromCell(r, 1:ncell(r))); colnames(xy) <- c("x","y")
  center_wi_1 <- c(2, -25)
  Sigma_wi_1  <- make_cov2D(angle_deg =  1, var_major = 2.0, var_minor = 2.0)
  omega_fun <- function(x, y){
    base_wi <- 0.04; cap_wi <- 0.90; amp1 <- 1.00
    val <- base_wi + amp1 * gauss2d(x, y, center_wi_1, Sigma_wi_1)
    pmax(0, pmin(val, cap_wi))
  }
  r[] <- with(xy, omega_fun(x, y)) * 0.4
  r
}

# Nombre de pixels par côté (grille spatiale pour la simulation)
n_pix_side <- 60

# Grilles d'effort d'échantillonnage
r_su <- build_omega_summer(xmin_all, xmax_all, ymin_all, ymax_all, nx = n_pix_side, ny = n_pix_side)
r_wi <- build_omega_winter(xmin_all, xmax_all, ymin_all, ymax_all, nx = n_pix_side, ny = n_pix_side)

# Effort d'échantillonnage associé à chaque individu
Omega_su_all <- terra::extract(r_su, X_su_all)[,1]
Omega_wi_all <- terra::extract(r_wi, X_wi_all)[,1]

# ======================= 3) Filtrage par Ω et table « long » =======================
# Chaque individu peut être observé en été ou en hiver
sel_ete   <- runif(N_tot) < Omega_su_all
sel_hiver <- runif(N_tot) < Omega_wi_all

obs_ete <- data.frame(
  id = which(sel_ete), saison = "ete",
  x_su = X_su_all[sel_ete,1], y_su = X_su_all[sel_ete,2],
  x_wi = X_wi_all[sel_ete,1], y_wi = X_wi_all[sel_ete,2],
  Omega = Omega_su_all[sel_ete]
)
obs_hiver <- data.frame(
  id = which(sel_hiver), saison = "hiver",
  x_su = X_su_all[sel_hiver,1], y_su = X_su_all[sel_hiver,2],
  x_wi = X_wi_all[sel_hiver,1], y_wi = X_wi_all[sel_hiver,2],
  Omega = Omega_wi_all[sel_hiver]
)

obs_long <- rbind(obs_ete, obs_hiver)
N <- nrow(obs_long)
X <- as.matrix(obs_long[,c("x_su","y_su","x_wi","y_wi")])
sampling.season <- ifelse(obs_long$saison=="ete", 1L, 2L)
Omega_vec <- obs_long$Omega

# ======================= 4) Génération des Omega dégradés =======================
# Fonction pour créer des Omega dégradés
create_degraded_omega <- function(r_original, fact_agg, noise_amp) {
  # Agrégation
  r_coarse <- terra::aggregate(r_original, fact = fact_agg, fun = mean, na.rm = TRUE)
  
  # Ajout du bruit
  vals <- values(r_coarse)[, 1]  
  n_cells <- length(vals)
  eps <- runif(n_cells, min = -noise_amp, max = noise_amp)
  vals_noisy <- vals * (1 + eps)
  vals_noisy <- pmax(0, pmin(1, vals_noisy))
  r_noisy <- r_coarse
  values(r_noisy) <- vals_noisy
  
  return(r_noisy)
}

# Créer les versions dégradées
r_su_deg <- create_degraded_omega(r_su, fact_agg = 4, noise_amp = 0.25)
r_wi_deg <- create_degraded_omega(r_wi, fact_agg = 4, noise_amp = 0.25)
r_su_deg2 <- create_degraded_omega(r_su, fact_agg = 8, noise_amp = 0.50)
r_wi_deg2 <- create_degraded_omega(r_wi, fact_agg = 8, noise_amp = 0.50)

# ======================= 5) Visualisation avec cowplot =======================
# Fonction pour convertir un SpatRaster en dataframe pour ggplot
raster_to_df <- function(r, name) {
  as.data.frame(r, xy = TRUE) %>% 
    setNames(c("x", "y", "value")) %>%
    mutate(raster_name = name)
}

# Convertir tous les rasters en dataframes
df_su <- raster_to_df(r_su, "Omega été original")
df_wi <- raster_to_df(r_wi, "Omega hiver original")
df_su_deg <- raster_to_df(r_su_deg, "Omega été dégradé")
df_wi_deg <- raster_to_df(r_wi_deg, "Omega hiver dégradé")
df_su_deg2 <- raster_to_df(r_su_deg2, "Omega été dégradé2")
df_wi_deg2 <- raster_to_df(r_wi_deg2, "Omega hiver dégradé2")

# Fonction pour créer un plot ggplot d'un raster
create_raster_plot <- function(df, title) {
  ggplot(df, aes(x = x, y = y, fill = value)) +
    geom_tile() +
    scale_fill_viridis_c(limits = c(0, 0.7), option = "plasma", name = "Omega") +
    labs(title = title, x = "X", y = "Y") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 10),
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 9),
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 8))
}

# Créer les 4 plots
p1 <- create_raster_plot(df_su, "Omega été original (60x60)")
p2 <- create_raster_plot(df_wi, "Omega hiver original (60x60)")
p3 <- create_raster_plot(df_su_deg, "Omega été dégradé (15x15)")
p4 <- create_raster_plot(df_wi_deg, "Omega hiver dégradé (15x15)")
p5 <- create_raster_plot(df_su_deg2, "Omega été dégradé2 (15x15)")
p6 <- create_raster_plot(df_wi_deg2, "Omega hiver dégradé2 (15x15)")

# Arranger les 4 plots en 2x2
plot_grid(p1, p2, p3, p4, p5, p6, ncol = 2, labels = "AUTO")

# ======================= 6) Ω dégradés au niveau des observations =======================
Omega_su_obs_deg <- terra::extract(r_su_deg, obs_ete[, c("x_su","y_su")])[,1]
Omega_wi_obs_deg <- terra::extract(r_wi_deg, obs_hiver[, c("x_wi","y_wi")])[,1]

Omega_vec_deg <- numeric(N)
Omega_vec_deg[obs_long$saison == "ete"]   <- Omega_su_obs_deg
Omega_vec_deg[obs_long$saison == "hiver"] <- Omega_wi_obs_deg