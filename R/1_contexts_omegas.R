# =============================================================
# Script 1/4 — Contextes biologiques (C1–C3) + Ω (original/deg1/deg2)
# =============================================================
# Ce script définit UNIQUEMENT les fonctions et objets nécessaires pour:
#  - décrire les 3 contextes biologiques (moyennes/covariances 4D)
#  - générer 3 versions d'effort d'échantillonnage Ω (été & hiver)
#  - quelques helpers de visualisation (optionnels)
#
# Les autres scripts devront simplement faire: source("1_contexts_omegas.R")
#
# Dépendances: Matrix, terra, ggplot2, dplyr, patchwork, viridisLite, stringr
# =============================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(terra)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(viridisLite)
  library(stringr)
})

# ---------- Utilitaires communs ----------
.gmmu_make_sigma <- function(sd, R){
  D4 <- diag(sd, 4, 4)
  S <- D4 %*% R %*% D4
  ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  if (!isSymmetric(S) || any(ev <= 1e-10)) S <- as.matrix(Matrix::nearPD(S)$mat)
  S
}

.gmmu_gauss2d <- function(x, y, mu, Sigma) {
  invS <- solve(Sigma)
  dx <- x - mu[1]; dy <- y - mu[2]
  Q <- invS[1,1]*dx*dx + (invS[1,2]+invS[2,1])*dx*dy + invS[2,2]*dy*dy
  exp(-0.5 * Q)
}

# ---------- Contextes biologiques (C1/C2/C3) ----------
# Retourne une liste: mu_true (list de K vecteurs 4D) et Sigma_true (list de K matrices 4x4)
get_context_params <- function(ctx = c("C1","C2","C3")) {
  ctx <- match.arg(ctx)
  if (ctx == "C1") {
    mu_true <- list(
      c(-4.0,  2.5, -4.0, -27.5),
      c( 0.0,  2.5,  0.0, -27.5),
      c( 4.0,  2.5,  4.0, -27.5)
    )
    R <- matrix(c(
      1.00,  0.10,  0.00,  0.00,
      0.10,  1.00,  0.00,  0.00,
      0.00,  0.00,  1.00,  0.10,
      0.00,  0.00,  0.10,  1.00
    ), 4, 4, byrow = TRUE)
    Sigma_true <- list(
      .gmmu_make_sigma(c(2,1,2,2), R),
      .gmmu_make_sigma(c(2,1,2,1), R),
      .gmmu_make_sigma(c(1,1,1,2), R)
    )
  }
  if (ctx == "C2") {
    mu_true <- list(
      c(-3.5,  0.0,  -1,   2),
      c( 0.0, -0.6,  -1,   2),
      c( 3.5,  0.0,   2,   5)
    )
    for (k in 1:3) mu_true[[k]][4] <- mu_true[[k]][4] - 30
    sd1 <- c(0.7, 2.8, 1.0, 2.0)
    R1 <- matrix(c(1,0,0,0, 0,1,0,0, 0,0,1,0.20, 0,0,0.20,1), 4, 4, byrow=TRUE)
    sd2 <- c(2.6, 2.4, 1.0, 2.0)
    R2 <- matrix(c(1,-0.70,0,0, -0.70,1,0,0, 0,0,1,0.20, 0,0,0.20,1), 4, 4, byrow=TRUE)
    sd3 <- c(0.8, 3.2, 0.6, 2.8)
    R3 <- diag(4)
    Sigma_true <- list(
      .gmmu_make_sigma(sd1, R1),
      .gmmu_make_sigma(sd2, R2),
      .gmmu_make_sigma(sd3, R3)
    )
  }
  if (ctx == "C3") {
    mu_true <- list(
      c(-3.5,  0.0,  0.2, -28.0),
      c( 0.0, -0.6,  0.3, -28.3),
      c( 3.5,  0.0,  0.0, -28.0)
    )
    sd1 <- c(0.7, 2.8, 0.7, 0.9)
    R1  <- matrix(c(1,0,0,0, 0,1,0,0, 0,0,1,0.10, 0,0,0.10,1), 4, 4, byrow=TRUE)
    sd2 <- c(2.6, 2.4, 0.7, 0.9)
    R2  <- matrix(c(1,-0.70,0,0, -0.70,1,0,0, 0,0,1,0.10, 0,0,0.10,1), 4, 4, byrow=TRUE)
    sd3 <- c(0.8, 3.2, 1.1, 1.6)
    R3  <- matrix(c(1,0,0,0, 0,1,0,0, 0,0,1,0.25, 0,0,0.25,1), 4, 4, byrow=TRUE)
    Sigma_true <- list(
      .gmmu_make_sigma(sd1, R1),
      .gmmu_make_sigma(sd2, R2),
      .gmmu_make_sigma(sd3, R3)
    )
  }
  list(mu_true = mu_true, Sigma_true = Sigma_true)
}

# ---------- Emprise commune (pour Ω et plots) ----------
.get_all_centers <- function(){
  lst <- lapply(c("C1","C2","C3"), function(ctx){
    pars <- get_context_params(ctx)
    do.call(rbind, pars$mu_true)
  })
  do.call(rbind, lst)
}
.centers_all <- .get_all_centers()
.EMP_PAD <- 10
.XMIN_ALL <- min(.centers_all[,c(1,3)]) - .EMP_PAD
.XMAX_ALL <- max(.centers_all[,c(1,3)]) + .EMP_PAD
.YMIN_ALL <- min(.centers_all[,c(2,4)]) - .EMP_PAD
.YMAX_ALL <- max(.centers_all[,c(2,4)]) + .EMP_PAD

# ---------- Ω (original + deux versions dégradées) ----------
build_omega_summer <- function(xmin=.XMIN_ALL, xmax=.XMAX_ALL, ymin=.YMIN_ALL, ymax=.YMAX_ALL, nx=60, ny=60){
  r <- rast(nrows = ny, ncols = nx, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
  xy <- as.data.frame(xyFromCell(r, 1:ncell(r))); colnames(xy) <- c("x","y")
  omega_fun <- function(x, y){
    base <- 0.04
    bumpL <- 0.90 * .gmmu_gauss2d(x, y, c(-3.8, 2.6), diag(c(0.60, 0.60)))
    bumpR <- 0.90 * .gmmu_gauss2d(x, y, c( 3.8, 2.6), diag(c(0.65, 0.65)))
    bumpC <- 1.00 * .gmmu_gauss2d(x, y, c( 0.0,-3.2), diag(c(6.0, 0.9)))
    pmax(0, pmin(base + bumpL + bumpR + bumpC, 0.70))
  }
  r[] <- with(xy, omega_fun(x, y)) * 0.5
  r
}

build_omega_winter <- function(xmin=.XMIN_ALL, xmax=.XMAX_ALL, ymin=.YMIN_ALL, ymax=.YMAX_ALL, nx=60, ny=60){
  r <- rast(nrows = ny, ncols = nx, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
  xy <- as.data.frame(xyFromCell(r, 1:ncell(r))); colnames(xy) <- c("x","y")
  center_wi_1 <- c(2, -25)
  Sigma_wi_1  <- diag(c(2.0, 2.0))
  omega_fun <- function(x, y){
    base_wi <- 0.04; cap_wi <- 0.90; amp1 <- 1.00
    val <- base_wi + amp1 * .gmmu_gauss2d(x, y, center_wi_1, Sigma_wi_1)
    pmax(0, pmin(val, cap_wi))
  }
  r[] <- with(xy, omega_fun(x, y)) * 0.4
  r
}

create_degraded_omega <- function(r_original, fact_agg, noise_amp) {
  r_coarse <- terra::aggregate(r_original, fact = fact_agg, fun = mean, na.rm = TRUE)
  vals <- values(r_coarse)[,1]
  eps <- runif(length(vals), min = -noise_amp, max = noise_amp)
  vals_noisy <- pmin(1, pmax(0, vals * (1 + eps)))
  r_noisy <- r_coarse
  values(r_noisy) <- vals_noisy
  r_noisy
}

# Génère une paire Ω été/hiver pour une version donnée
make_omega_pair <- function(version=c("original","deg1","deg2"), n_pix_side=60){
  version <- match.arg(version)
  r_su <- build_omega_summer(nx=n_pix_side, ny=n_pix_side)
  r_wi <- build_omega_winter(nx=n_pix_side, ny=n_pix_side)
  if (version=="original") return(list(summer=r_su, winter=r_wi))
  if (version=="deg1") return(list(summer=create_degraded_omega(r_su, 4, 0.25),
                                   winter=create_degraded_omega(r_wi, 4, 0.25)))
  if (version=="deg2") return(list(summer=create_degraded_omega(r_su, 8, 0.50),
                                   winter=create_degraded_omega(r_wi, 8, 0.50)))
}

# ---------- Helpers d'inspection (facultatifs) ----------
plot_omega_pair <- function(omega_pair, title_prefix="Ω"){
  df_from_r <- function(r, title){ as.data.frame(r, xy=TRUE) |>
      setNames(c("x","y","value")) |>
      mutate(title=title) }
  p <- function(df){ ggplot(df, aes(x,y,fill=value)) + geom_tile() +
      scale_fill_viridis_c(limits=c(0,0.7), option="plasma", name="Omega") +
      coord_fixed(xlim=c(.XMIN_ALL,.XMAX_ALL), ylim=c(.YMIN_ALL,.YMAX_ALL), expand=FALSE) +
      theme_minimal() + labs(title=unique(df$title), x="X", y="Y") }
  (p(df_from_r(omega_pair$summer, sprintf("%s été", title_prefix))) /
      p(df_from_r(omega_pair$winter, sprintf("%s hiver", title_prefix))))
}

# ---------- Plot final de vérification ----------
# omega_pair_demo <- make_omega_pair("deg1", n_pix_side=60)
# print(plot_omega_pair(omega_pair_demo, "Ω deg1"))
# 
# # Vérification visuelle des ellipses du contexte C2 (par ex.)
# ctx_demo <- get_context_params("C2")
# print(plot_ctx_truth("C2"))


# ---- FIN Script 1 ----
