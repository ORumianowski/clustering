# ==================================================
# Simulation 4D (mu & Sigma "en une fois") + Ω été (original) + Ω hiver (explicite large)
# + filtrage + plots de contrôle (modifs: hiver y-10 ; emprise Ω commune ; plots combinés)
# ==================================================

rm(list = ls())
set.seed(123)

suppressPackageStartupMessages({
  library(MASS)
  library(Matrix)   # nearPD
  library(terra)
  library(ggplot2)
  library(patchwork)
})

# -------------------- utilitaires --------------------
gauss2d <- function(x, y, mu, Sigma) {
  invS <- solve(Sigma)
  dx <- x - mu[1]; dy <- y - mu[2]
  Q <- invS[1,1]*dx*dx + (invS[1,2]+invS[2,1])*dx*dy + invS[2,2]*dy*dy
  exp(-0.5 * Q)
}
build_sigma <- function(sd, R){
  D <- diag(sd, 4, 4)
  S <- D %*% R %*% D
  ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  if (!isSymmetric(S) || any(ev <= 1e-10)) S <- as.matrix(nearPD(S)$mat)
  S
}
make_cov2D <- function(angle_deg, var_major, var_minor){
  phi <- angle_deg * pi/180
  R <- matrix(c(cos(phi), -sin(phi), sin(phi), cos(phi)), 2, 2, byrow = TRUE)
  S <- R %*% diag(c(var_major, var_minor)) %*% t(R)
  ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  if (!isSymmetric(S) || any(ev <= 1e-10)) S <- as.matrix(nearPD(S)$mat)
  S
}

# -------------------- paramètres GMM 4D --------------------
K <- 3
N_tot <- 5000
pi_true <- c(0.20, 0.45, 0.35)

# Moyennes 4D (lon_été, lat_été, lon_hiver, lat_hiver) — stéréotypées
mu <- list(
  c(-3.5,  0.0,  -3.0,  1.2),  # C1 : été = gauche, hiver = gauche (diag +45°)
  c( 0.0, -0.6,   0.2, -0.4),  # C2 : été = diagonal -45° (centre), hiver = ~vertical centre
  c( 3.5,  0.0,   3.2,  0.3)   # C3 : été = droite (vertical), hiver = droite (horizontal)
)

# (i) DÉCALAGE GLOBAL HIVER : y_hiver - 10
for(k in 1:K) mu[[k]][4] <- mu[[k]][4] - 30

# Covariances 4x4 complètes (corrélations été↔hiver) — stéréotypées
## C1 : Été vertical (sd_x petit, sd_y grand), Hiver diagonal +45°
sd1 <- c(0.7, 2.8, 2.6, 1.9)  # (x_su, y_su, x_wi, y_wi)
R1 <- matrix(c(
  1.00,  0.00,  0.55, -0.15,
  0.00,  1.00,  0.10,  0.30,
  0.55,  0.10,  1.00,  0.65,
  -0.15,  0.30,  0.65,  1.00
), 4, 4, byrow = TRUE)
Sigma1 <- build_sigma(sd1, R1)

## C2 : Été diagonal -45° (corr < 0), Hiver vertical (sd_x petit, sd_y grand)
sd2 <- c(2.6, 2.4, 0.8, 3.0)
R2 <- matrix(c(
  1.00, -0.70,  0.40, -0.35,
  -0.70,  1.00, -0.10,  0.45,
  0.40, -0.10,  1.00,  0.00,
  -0.35,  0.45,  0.00,  1.00
), 4, 4, byrow = TRUE)
Sigma2 <- build_sigma(sd2, R2)

## C3 : Été vertical (à droite), Hiver horizontal (sd_x grand, sd_y petit)
sd3 <- c(0.8, 3.2, 3.1, 0.9)
R3 <- matrix(c(
  1.00,  0.00,  0.50,  0.05,
  0.00,  1.00, -0.20,  0.25,
  0.50, -0.20,  1.00,  0.05,
  0.05,  0.25,  0.05,  1.00
), 4, 4, byrow = TRUE)
Sigma3 <- build_sigma(sd3, R3)

Sigma <- list(Sigma1, Sigma2, Sigma3)

# -------------------- simulation 4D (sans Ω) --------------------
z <- sample(1:K, N_tot, replace = TRUE, prob = pi_true)
X4_all <- t(vapply(
  1:N_tot, function(i) MASS::mvrnorm(1, mu = mu[[z[i]]], Sigma = Sigma[[z[i]]]),
  numeric(4)
))
colnames(X4_all) <- c("x_su","y_su","x_wi","y_wi")

# Séparations 2D pour Ω/plots
X_su_all <- X4_all[, 1:2, drop = FALSE]
X_wi_all <- X4_all[, 3:4, drop = FALSE]

# -------------------- Ω(x, saison) --------------------
# ÉTÉ : identique à l’original
build_omega_summer <- function(xmin, xmax, ymin, ymax, nx = 60, ny = 60){
  r <- rast(nrows = ny, ncols = nx, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
  xy <- as.data.frame(xyFromCell(r, 1:ncell(r))); colnames(xy) <- c("x","y")
  omega_fun <- function(x, y){
    base <- 0.04
    bumpL <- 0.90 * gauss2d(x, y, c(-3.8, 2.6), diag(c(0.60, 0.60)))
    bumpR <- 0.90 * gauss2d(x, y, c( 3.8, 2.6), diag(c(0.65, 0.65)))
    bumpC <- 1.00 * gauss2d(x, y, c( 0.0,-3.2), diag(c(6.0, 0.9)))  # “trou” bas large
    pmax(0, pmin(base + bumpL + bumpR + bumpC, 0.70))
  }
  r[] <- with(xy, omega_fun(x, y))
  r
}

# HIVER : 2 larges taches déterministes — centrées près des clusters d’hiver
# (ii) Comme les centres dépendent de mu (déjà y-10), le décalage est automatique.
build_omega_winter_explicit <- function(xmin, xmax, ymin, ymax, nx = 80, ny = 80,
                                        mu, enlarge = 1.0){
  center_wi_1 <- c(mu[[1]][3], mu[[1]][4]) + c( 0.6, -1.5)
  # center_wi_2 <- c(mu[[3]][3], mu[[3]][4]) + c(-2.5,  -5.4)
  Sigma_wi_1 <- make_cov2D(angle_deg =  01, var_major = 5.0*enlarge, var_minor = 2.0*enlarge)
  # Sigma_wi_2 <- make_cov2D(angle_deg =  -5, var_major =  5.0*enlarge, var_minor = 1.6*enlarge)
  base_wi <- 0.04; cap_wi <- 0.90; amps_wi <- c(1.00, 0.95)
  
  r <- rast(nrows = ny, ncols = nx, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
  xy <- as.data.frame(xyFromCell(r, 1:ncell(r))); colnames(xy) <- c("x","y")
  omega_fun <- function(x, y){
    val <- base_wi
    val <- val + amps_wi[1] * gauss2d(x, y, center_wi_1, Sigma_wi_1)
    # val <- val + amps_wi[2] * gauss2d(x, y, center_wi_2, Sigma_wi_2)
    pmax(0, pmin(val, cap_wi))
  }
  r[] <- with(xy, omega_fun(x, y))
  r
}

# -------------------- Emprise commune + rasters Ω --------------------
pad <- 1
xmin_all <- min(c(X_su_all[,1], X_wi_all[,1])) - pad
xmax_all <- max(c(X_su_all[,1], X_wi_all[,1])) + pad
ymin_all <- min(c(X_su_all[,2], X_wi_all[,2])) - pad
ymax_all <- max(c(X_su_all[,2], X_wi_all[,2])) + pad

# (ii) même extension spatiale pour Ω été et Ω hiver
r_su <- build_omega_summer(xmin_all, xmax_all, ymin_all, ymax_all, nx = 80, ny = 80)
r_wi <- build_omega_winter_explicit(xmin_all, xmax_all, ymin_all, ymax_all, nx = 80, ny = 80, mu = mu)

# -------------------- filtrage par Ω --------------------
Omega_su_all <- terra::extract(r_su, X_su_all, method = "simple")[,1]
Omega_wi_all <- terra::extract(r_wi, X_wi_all, method = "simple")[,1]
keep_su <- runif(N_tot) < Omega_su_all
keep_wi <- runif(N_tot) < Omega_wi_all

X_su_obs <- X_su_all[keep_su, , drop = FALSE]
X_wi_obs <- X_wi_all[keep_wi, , drop = FALSE]

cat("N observés ÉTÉ  après Ω :", nrow(X_su_obs), "\n")
cat("N observés HIVER après Ω :", nrow(X_wi_obs), "\n")

# -------------------- plots de contrôle --------------------

theme0 <- theme_minimal(base_size = 12)
square_theme <- theme0 + theme(aspect.ratio = 1)  # évite tout étirement

# AVANT Ω — points combinés (été = orange, hiver = bleu foncé)
before_df <- rbind(
  data.frame(x = X_su_all[,1], y = X_su_all[,2], saison = "été"),
  data.frame(x = X_wi_all[,1], y = X_wi_all[,2], saison = "hiver")
)
p_before <- ggplot(before_df, aes(x = x, y = y, color = saison)) +
  geom_point(size = 0.55, alpha = 0.85) +
  scale_color_manual(values = c("été" = "orange", "hiver" = "darkblue")) +
  coord_equal(xlim = c(xmin_all, xmax_all), ylim = c(ymin_all, ymax_all), expand = FALSE) +
  square_theme + labs(title = "AVANT Ω — Été & Hiver", x = "x", y = "y", color = "Saison")

# Ω ÉTÉ (emprise commune)
su_r_df <- as.data.frame(r_su, xy = TRUE); names(su_r_df)[3] <- "omega"
p_su_omega <- ggplot(su_r_df, aes(x = x, y = y, fill = omega)) +
  geom_raster() +
  coord_equal(xlim = c(xmin_all, xmax_all), ylim = c(ymin_all, ymax_all), expand = FALSE) +
  scale_fill_viridis_c(name = "Ω (été)") +
  square_theme + labs(title = "Fond Ω — ÉTÉ", x = "x", y = "y")

# Ω HIVER (emprise commune)
wi_r_df <- as.data.frame(r_wi, xy = TRUE); names(wi_r_df)[3] <- "omega"
p_wi_omega <- ggplot(wi_r_df, aes(x = x, y = y, fill = omega)) +
  geom_raster() +
  coord_equal(xlim = c(xmin_all, xmax_all), ylim = c(ymin_all, ymax_all), expand = FALSE) +
  scale_fill_viridis_c(name = "Ω (hiver)") +
  square_theme + labs(title = "Fond Ω — HIVER", x = "x", y = "y")

# APRÈS Ω — points combinés
after_df <- rbind(
  data.frame(x = X_su_obs[,1], y = X_su_obs[,2], saison = "été"),
  data.frame(x = X_wi_obs[,1], y = X_wi_obs[,2], saison = "hiver")
)
p_after <- ggplot(after_df, aes(x = x, y = y, color = saison)) +
  geom_point(size = 0.55, alpha = 0.9) +
  scale_color_manual(values = c("été" = "orange", "hiver" = "darkblue")) +
  coord_equal(xlim = c(xmin_all, xmax_all), ylim = c(ymin_all, ymax_all), expand = FALSE) +
  square_theme + labs(title = "APRÈS Ω — Été & Hiver", x = "x", y = "y", color = "Saison")

# -------------------- Disposition finale demandée --------------------
(p_before | p_after) / (p_su_omega | p_wi_omega)
