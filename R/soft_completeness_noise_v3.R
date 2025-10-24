# ==================================================
# Contextes biologiques (ellipses soignées) + Ω (3 versions)
# AUCUNE application d'Ω aux contextes (juste des figures)
# ==================================================

rm(list = ls())
set.seed(123)

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(Matrix)
  library(terra); library(viridisLite); library(patchwork); library(stringr)
})

# ----------------------- Utilitaires généraux -----------------------
build_sigma <- function(sd, R){
  D4 <- diag(sd, 4, 4)
  S <- D4 %*% R %*% D4
  ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  if (!isSymmetric(S) || any(ev <= 1e-10)) S <- as.matrix(Matrix::nearPD(S)$mat)
  S
}

# Ellipse lissée (niveaux 50/75/90 %)
ellipse_points <- function(mu, Sigma, r, n = 300){
  ang <- seq(0, 2*pi, length.out = n)
  circle <- cbind(cos(ang), sin(ang))
  R <- chol(Sigma)
  pts <- circle %*% t(R) * r
  sweep(pts, 2, mu, "+")
}
make_aplats_sets <- function(mu_su, Sig_su, mu_wi, Sig_wi,
                             probs = c(0.50, 0.75, 0.90)){
  r_levels <- sqrt(qchisq(probs, df = 2))
  cov_labels <- sprintf("%d%%", round(100*probs))
  out <- vector("list", 2*nrow(mu_su)*length(r_levels)); idx <- 1L
  for (k in 1:nrow(mu_su)){
    for (i in seq_along(r_levels)){
      pts <- ellipse_points(mu_su[k,], Sig_su[[k]], r_levels[i])
      out[[idx]] <- data.frame(x=pts[,1], y=pts[,2],
                               cluster=factor(k), set="été", level=cov_labels[i]); idx <- idx+1L
    }
    for (i in seq_along(r_levels)){
      pts <- ellipse_points(mu_wi[k,], Sig_wi[[k]], r_levels[i])
      out[[idx]] <- data.frame(x=pts[,1], y=pts[,2],
                               cluster=factor(k), set="hiver", level=cov_labels[i]); idx <- idx+1L
    }
  }
  out <- do.call(rbind, out)
  out$set   <- factor(out$set, levels=c("été","hiver"))
  out$level <- factor(out$level, levels=c("50%","75%","90%"))
  out
}

# ----------------------- Contextes biologiques -----------------------
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
      build_sigma(c(2,1,2,2), R),
      build_sigma(c(2,1,2,1), R),
      build_sigma(c(1,1,1,2), R)
    )
  }
  
  if (ctx == "C2") {  # == ton paramétrage actuel
    mu_true <- list(
      c(-3.5,  0.0,  -1,   2),
      c( 0.0, -0.6,  -1,   2),
      c( 3.5,  0.0,   2,   5)
    )
    for (k in 1:3) mu_true[[k]][4] <- mu_true[[k]][4] - 30
    
    sd1 <- c(0.7, 2.8, 1.0, 2.0)
    R1 <- matrix(c(
      1.00,  0.00,  0.00,  0.00,
      0.00,  1.00,  0.00,  0.00,
      0.00,  0.00,  1.00,  0.20,
      0.00,  0.00,  0.20,  1.00
    ), 4, 4, byrow = TRUE)
    
    sd2 <- c(2.6, 2.4, 1.0, 2.0)
    R2 <- matrix(c(
      1.00, -0.70,  0.00,  0.00,
      -0.70,  1.00,  0.00,  0.00,
      0.00,  0.00,  1.00,  0.20,
      0.00,  0.00,  0.20,  1.00
    ), 4, 4, byrow = TRUE)
    
    sd3 <- c(0.8, 3.2, 0.6, 2.8)
    R3 <- diag(4)
    
    Sigma_true <- list(
      build_sigma(sd1, R1),
      build_sigma(sd2, R2),
      build_sigma(sd3, R3)
    )
  }
  
  if (ctx == "C3") {
    mu_true <- list(
      c(-2.2,  1.0,  0.2, -28.0),
      c(-0.8,  0.6,  0.3, -28.3),
      c( 3.2,  1.2,  0.0, -28.0)
    )
    R_k12 <- matrix(c(
      1.00, -0.75,  0.00,  0.00,
      -0.75, 1.00,  0.00,  0.00,
      0.00,  0.00,  1.00,  0.10,
      0.00,  0.00,  0.10,  1.00
    ), 4, 4, byrow = TRUE)
    sd_k12 <- c(2.4, 2.2, 0.7, 0.9)
    
    R_k3 <- matrix(c(
      1.00,  0.15,  0.00,  0.00,
      0.15,  1.00,  0.00,  0.00,
      0.00,  0.00,  1.00,  0.25,
      0.00,  0.00,  0.25,  1.00
    ), 4, 4, byrow = TRUE)
    sd_k3 <- c(2, 1.3, 1.1, 1.6)
    
    Sigma_true <- list(
      build_sigma(sd_k12, R_k12),
      build_sigma(sd_k12, R_k12),
      build_sigma(sd_k3,  R_k3 )
    )
  }
  
  list(mu_true = mu_true, Sigma_true = Sigma_true)
}

# ----------------------- Plot « Script 2 » appliqué à la vérité only -----------------------
plot_ctx_truth <- function(ctx = c("C1","C2","C3"), panel_aspect = 0.9){
  ctx <- match.arg(ctx)
  pars <- get_context_params(ctx)
  mu_true    <- pars$mu_true
  Sigma_true <- pars$Sigma_true
  K <- length(mu_true)
  
  mu_mat <- do.call(rbind, mu_true)
  mu_su  <- mu_mat[,1:2,drop=FALSE]
  mu_wi  <- mu_mat[,3:4,drop=FALSE]
  Sig_su <- lapply(Sigma_true, \(S) S[1:2,1:2,drop=FALSE])
  Sig_wi <- lapply(Sigma_true, \(S) S[3:4,3:4,drop=FALSE])
  
  polys_true <- make_aplats_sets(mu_su, Sig_su, mu_wi, Sig_wi)
  
  # segments été→hiver: on relie les centres (mu)
  segs <- data.frame(
    x = mu_su[,1], y = mu_su[,2],
    xend = mu_wi[,1], yend = mu_wi[,2],
    cluster = factor(1:K)
  )
  
  # emprise jolie
  pad <- 1
  xlim_all <- range(c(mu_su[,1], mu_wi[,1])) + c(-pad, pad)
  ylim_all <- range(c(mu_su[,2], mu_wi[,2])) + c(-pad, pad)
  
  # palettes & alphas
  cluster_palette <- setNames(rocket(K, begin=0.10, end=0.75, direction=-1), as.character(1:K))
  coverage_alphas <- c("50%"=.50, "75%"=.20, "90%"=.10)
  
  ggplot() +
    geom_segment(data = segs,
                 aes(x=x, y=y, xend=xend, yend=yend, color=cluster),
                 alpha=.8, linewidth=1) +
    geom_polygon(data = polys_true,
                 aes(x=x, y=y, group=interaction(cluster, level, set),
                     fill=cluster, alpha=level), color=NA) +
    scale_fill_manual(values=cluster_palette, name="Clusters",
                      guide = guide_legend(override.aes = list(alpha = 1))) +
    scale_color_manual(values=cluster_palette, guide = "none") +
    scale_alpha_manual(values=coverage_alphas, name="Contours prob.") +
    coord_fixed(ratio = panel_aspect, xlim=xlim_all, ylim=ylim_all, expand=FALSE) +
    theme_minimal(base_size=12) +
    theme(plot.title = element_text(face="bold"),
          legend.position = "bottom") +
    labs(title = paste0("Ellipses combinées — VRAI (", ctx, ")"),
         subtitle = "Segments = centres été → hiver",
         x="x", y="y")
}

# ----------------------- Ω : original + 2 dégradés (été & hiver) -----------------------
gauss2d <- function(x, y, mu, Sigma) {
  invS <- solve(Sigma)
  dx <- x - mu[1]; dy <- y - mu[2]
  Q <- invS[1,1]*dx*dx + (invS[1,2]+invS[2,1])*dx*dy + invS[2,2]*dy*dy
  exp(-0.5 * Q)
}
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
  Sigma_wi_1  <- diag(c(2.0, 2.0))
  omega_fun <- function(x, y){
    base_wi <- 0.04; cap_wi <- 0.90; amp1 <- 1.00
    val <- base_wi + amp1 * gauss2d(x, y, center_wi_1, Sigma_wi_1)
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
raster_to_df <- function(r, title) {
  as.data.frame(r, xy = TRUE) %>%
    setNames(c("x","y","value")) %>%
    mutate(title = title)
}
plot_raster <- function(df) {
  ggplot(df, aes(x, y, fill = value)) +
    geom_tile() +
    scale_fill_viridis_c(limits = c(0, 0.7), option = "plasma", name = "Omega") +
    labs(title = unique(df$title), x = "X", y = "Y") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 11))
}

# ----------------------- FIGURES -----------------------
# 1) Ellipses jolies pour les 3 contextes
pC1 <- plot_ctx_truth("C1")
pC2 <- plot_ctx_truth("C2")  # ton cas « milieu »
pC3 <- plot_ctx_truth("C3")
(pC1 | pC2 | pC3) + plot_layout(heights = c(1,1,1))

# 2) Ω : original + dégradé1 + dégradé2 (été & hiver)
#    — emprise large couvrant tous les centres
get_all_centers <- function(){
  lapply(c("C1","C2","C3"), function(ctx){
    pars <- get_context_params(ctx)
    do.call(rbind, pars$mu_true)
  }) %>% do.call(rbind, .)
}
centers_all <- get_all_centers()
pad <- 2
xmin_all <- min(centers_all[,c(1,3)]) - pad
xmax_all <- max(centers_all[,c(1,3)]) + pad
ymin_all <- min(centers_all[,c(2,4)]) - pad
ymax_all <- max(centers_all[,c(2,4)]) + pad

n_pix_side <- 60
r_su      <- build_omega_summer(xmin_all, xmax_all, ymin_all, ymax_all, nx = n_pix_side, ny = n_pix_side)
r_wi      <- build_omega_winter(xmin_all, xmax_all, ymin_all, ymax_all, nx = n_pix_side, ny = n_pix_side)
r_su_deg1 <- create_degraded_omega(r_su, fact_agg = 4, noise_amp = 0.25)
r_wi_deg1 <- create_degraded_omega(r_wi, fact_agg = 4, noise_amp = 0.25)
r_su_deg2 <- create_degraded_omega(r_su, fact_agg = 8, noise_amp = 0.50)
r_wi_deg2 <- create_degraded_omega(r_wi, fact_agg = 8, noise_amp = 0.50)

df_su      <- raster_to_df(r_su,      "Ω été — original")
df_wi      <- raster_to_df(r_wi,      "Ω hiver — original")
df_su_d1   <- raster_to_df(r_su_deg1, "Ω été — dégradé 1 (x4, ±25%)")
df_wi_d1   <- raster_to_df(r_wi_deg1, "Ω hiver — dégradé 1 (x4, ±25%)")
df_su_d2   <- raster_to_df(r_su_deg2, "Ω été — dégradé 2 (x8, ±50%)")
df_wi_d2   <- raster_to_df(r_wi_deg2, "Ω hiver — dégradé 2 (x8, ±50%)")

p1 <- plot_raster(df_su)
p2 <- plot_raster(df_wi)
p3 <- plot_raster(df_su_d1)
p4 <- plot_raster(df_wi_d1)
p5 <- plot_raster(df_su_d2)
p6 <- plot_raster(df_wi_d2)

(p1 | p2) / (p3 | p4) / (p5 | p6)
