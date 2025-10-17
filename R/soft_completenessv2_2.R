# ==================================================
# GMM Bayésien 4D (NIMBLE) avec correction Ω_saisons :
#   log p_o(x|θ,Ω) = log Ω_su(x_su) + log Ω_wi(x_wi) + log p(x|θ) - log Z(θ,Ω)
# Z(θ,Ω) évalué par MC sur une grille produit (été × hiver)
#   ones-trick : 1 ~ Bernoulli(exp(ℓ_i - C_ub))
# ==================================================

rm(list = ls())
set.seed(123)

suppressPackageStartupMessages({
  library(MASS)
  library(Matrix)   # nearPD
  library(terra)
  library(ggplot2)
  library(patchwork)
  library(nimble)
  library(mvtnorm)
})

# -------------------- utilitaires --------------------
gauss2d <- function(x, y, mu, Sigma) {
  invS <- solve(Sigma)
  dx <- x - mu[1]; dy <- y - mu[2]
  Q <- invS[1,1]*dx*dx + (invS[1,2]+invS[2,1])*dx*dy + invS[2,2]*dy*dy
  exp(-0.5 * Q)
}
build_sigma <- function(sd, R){
  D4 <- diag(sd, 4, 4)
  S <- D4 %*% R %*% D4
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
D <- 4
N_tot <- 5000
pi_true <- c(0.20, 0.45, 0.35)

# Moyennes 4D (lon_été, lat_été, lon_hiver, lat_hiver) — stéréotypées
mu <- list(
  c(-3.5,  0.0,  -3.0,  1.2),  # C1 : été = gauche (vertical), hiver = gauche (diag +45°)
  c( 0.0, -0.6,   0.2, -0.4),  # C2 : été = diagonal -45° (centre), hiver = ~vertical centre
  c( 3.5,  0.0,   3.2,  0.3)   # C3 : été = droite (vertical), hiver = droite (horizontal)
)
# Décalage global hiver (axe y) pour séparer visuellement les saisons
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
true_Sigma <- Sigma  # pour les plots

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
# ÉTÉ : identique à l’original (3 bosses)
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
  r[] <- with(xy, omega_fun(x, y))
  r
}

# HIVER : 1 large tache explicite proche du cluster hiver C1 (décalée pour limiter le recouvrement)
build_omega_winter_explicit <- function(xmin, xmax, ymin, ymax, nx = 80, ny = 80,
                                        mu, enlarge = 1.0){
  center_wi_1 <- c(mu[[1]][3], mu[[1]][4]) + c( 0.6, -1.5)
  Sigma_wi_1  <- make_cov2D(angle_deg =  1, var_major = 5.0*enlarge, var_minor = 2.0*enlarge)
  base_wi <- 0.04; cap_wi <- 0.90; amp1 <- 1.00
  
  r <- rast(nrows = ny, ncols = nx, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
  xy <- as.data.frame(xyFromCell(r, 1:ncell(r))); colnames(xy) <- c("x","y")
  omega_fun <- function(x, y){
    val <- base_wi
    val <- val + amp1 * gauss2d(x, y, center_wi_1, Sigma_wi_1)
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

# même extension spatiale pour Ω été et Ω hiver
r_su <- build_omega_summer(xmin_all, xmax_all, ymin_all, ymax_all, nx = 80, ny = 80)
r_wi <- build_omega_winter_explicit(xmin_all, xmax_all, ymin_all, ymax_all, nx = 80, ny = 80, mu = mu)

# -------------------- filtrage par Ω_total = Ω_été * Ω_hiver --------------------
Omega_su_all <- terra::extract(r_su, X_su_all, method = "simple")[,1]
Omega_wi_all <- terra::extract(r_wi, X_wi_all, method = "simple")[,1]
Omega_all    <- pmax(0, pmin(1, Omega_su_all * Omega_wi_all))

keep <- runif(N_tot) < Omega_all
X <- X4_all[keep,,drop=FALSE]; N <- nrow(X)
Omega_vec <- Omega_all[keep]
cat("N observés après Ω_total :", N, "\n")

# ----------------------- 3) Grilles pour Z(θ,Ω) (MC) -----------------------
# On utilise deux grilles 2D (été & hiver) et on échantillonne Mmc paires (produit)
su_df <- as.data.frame(r_su, xy=TRUE); names(su_df)[3] <- "omega_su"
wi_df <- as.data.frame(r_wi, xy=TRUE); names(wi_df)[3] <- "omega_wi"
A_su <- prod(res(r_su)); A_wi <- prod(res(r_wi))  # aires de cellule 2D
A_tot <- A_su * A_wi

set.seed(1234)
Mmc <- 6000L
idx_su <- sample.int(nrow(su_df), Mmc, replace=TRUE)
idx_wi <- sample.int(nrow(wi_df), Mmc, replace=TRUE)

grid_su <- as.matrix(su_df[idx_su, c("x","y")])
grid_wi <- as.matrix(wi_df[idx_wi, c("x","y")])
omega_su_mc <- su_df$omega_su[idx_su]
omega_wi_mc <- wi_df$omega_wi[idx_wi]

# ----------------------- 4) Modèle corrigé Ω (4D) -----------------------
dmvnorm_nimble <- nimbleFunction(
  run = function(x = double(1), mean = double(1), Prec = double(2),
                 log = logical(0, default = FALSE)) {
    returnType(double(0))
    Dloc <- length(x); xm <- x - mean
    qf <- inprod(xm, Prec %*% xm)
    U <- chol(Prec)
    ldet <- 2 * sum(log(diag(U)))
    logdens <- 0.5 * ldet - 0.5 * Dloc * log(2*pi) - 0.5 * qf
    if (log) return(logdens) else return(exp(logdens))
  }
)

# log Z(θ,Ω) par MC sur couples (grid_su[m,], grid_wi[m,])
logZ_calc_mc <- nimbleFunction(
  run = function(mu   = double(2),   # K x D
                 Prec = double(3),   # D x D x K
                 pi   = double(1),   # K
                 grid_su = double(2),# Mmc x 2
                 grid_wi = double(2),# Mmc x 2
                 omega_su = double(1), # Mmc
                 omega_wi = double(1), # Mmc
                 A_tot = double(0),
                 K = integer(0), D = integer(0), Mmc = integer(0)) {
    returnType(double(0))
    sumZ <- 0.0
    for(m in 1:Mmc){
      mix <- 0.0
      for(k in 1:K){
        x4 <- c(grid_su[m,1], grid_su[m,2], grid_wi[m,1], grid_wi[m,2])
        mix <- mix + pi[k] * dmvnorm_nimble(x4, mu[k,1:D], Prec[1:D,1:D,k], FALSE)
      }
      sumZ <- sumZ + omega_su[m] * omega_wi[m] * mix
    }
    Z <- A_tot * (sumZ / Mmc)
    if (Z < 1e-300) Z <- 1e-300
    return(log(Z))
  })

C_ub   <- 200
p_min  <- 1e-300
p_max1 <- 1 - 1e-12

code_corrige <- nimbleCode({
  logZ <- logZ_calc_mc(mu[1:K,1:D], Prec[1:D,1:D,1:K], pi[1:K],
                       grid_su[1:Mmc,1:2], grid_wi[1:Mmc,1:2],
                       omega_su[1:Mmc], omega_wi[1:Mmc],
                       A_tot, K, D, Mmc)
  for (i in 1:N) {
    for (k in 1:K) {
      dens[i,k] <- pi[k] * exp(dmvnorm_nimble(X[i,1:D], mu[k,1:D], Prec[1:D,1:D,k], TRUE))
    }
    mixdens[i] <- sum(dens[i,1:K])
    ll[i] <- log(Omega[i] * mixdens[i] + 1e-300) - logZ
    p_raw[i]      <- exp(ll[i] - C_ub)
    p_clip_hi[i]  <- min(p_raw[i], p_max1)
    p[i]          <- max(p_clip_hi[i], p_min)
    ones[i] ~ dbern(p[i])
  }
  pi[1:K] ~ ddirch(alpha[1:K])
  for (k in 1:K) {
    mu[k,1:D] ~ dmnorm(mu0[1:D], prec = Prec0[1:D,1:D])
    Prec[1:D,1:D,k] ~ dwish(R[1:D,1:D], df)
  }
})

# constants / data
Prec0 <- diag(D) * 1e-2
R     <- diag(D)
constants_corr <- list(
  N=N, D=D, K=K,
  Mmc=Mmc, A_tot=A_tot,
  alpha=rep(1,K),
  mu0=rep(0,D),
  Prec0=Prec0,
  R=R, df=D+2,
  p_min=p_min, p_max1=p_max1, C_ub=C_ub
)
data_corr <- list(
  X = X,
  Omega = Omega_vec,
  ones = rep(1, N),
  grid_su = grid_su,
  grid_wi = grid_wi,
  omega_su = omega_su_mc,
  omega_wi = omega_wi_mc
)

# inits
Prec_init <- array(0, dim=c(D,D,K)); for(k in 1:K) Prec_init[,,k] <- diag(D)
inits_corr <- list(mu = matrix(rnorm(K*D), K, D), Prec = Prec_init, pi = rep(1/K,K))

# build & run
model_corr  <- nimbleModel(code_corrige, data=data_corr, constants=constants_corr, inits=inits_corr, check=FALSE)
cmodel_corr <- compileNimble(model_corr)
conf_corr   <- configureMCMC(model_corr, monitors=c("mu","pi","Prec"))
mcmc_corr   <- buildMCMC(conf_corr)
cmcmc_corr  <- compileNimble(mcmc_corr, project=model_corr)
samples_corr <- runMCMC(cmcmc_corr, niter=3000, nburnin=1500, thin=6, setSeed=TRUE)

# ----------------------- 5) Modèle naïf (sans Ω) en 4D -----------------------
code_naif <- nimbleCode({
  for (i in 1:N){
    z[i] ~ dcat(pi[1:K])
    X[i,1:D] ~ dmnorm(mu[z[i],1:D], prec=Prec[1:D,1:D,z[i]])
  }
  pi[1:K] ~ ddirch(alpha[1:K])
  for (k in 1:K){
    mu[k,1:D] ~ dmnorm(mu0[1:D], prec=Prec0[1:D,1:D])
    Prec[1:D,1:D,k] ~ dwish(R[1:D,1:D], df)
  }
})
constants_naif <- list(N=N,D=D,K=K,alpha=rep(1,K),mu0=rep(0,D),Prec0=Prec0,R=R,df=D+2)
data_naif <- list(X=X)
inits_naif <- list(z=sample(1:K,N,TRUE), mu=matrix(rnorm(K*D),K,D), Prec=Prec_init, pi=rep(1/K,K))
model_naif  <- nimbleModel(code_naif, data=data_naif, constants=constants_naif, inits=inits_naif, check=FALSE)
cmodel_naif <- compileNimble(model_naif)
conf_naif   <- configureMCMC(model_naif, monitors=c("mu","pi","Prec"))
mcmc_naif   <- buildMCMC(conf_naif)
cmcmc_naif  <- compileNimble(mcmc_naif, project=model_naif)
samples_naif <- runMCMC(cmcmc_naif, niter=4000, nburnin=2000, thin=8, setSeed=TRUE)

# ----------------------- 6) Plots (4D -> 2D été/hiver) -----------------------
# Helpers
extr_muSigma4D <- function(samples, K, D = 4){
  mu_cols <- grep("^mu\\[", colnames(samples), value = TRUE)
  k_idx <- as.integer(sub("^mu\\[(\\d+),\\s*(\\d+)\\]$", "\\1", mu_cols))
  d_idx <- as.integer(sub("^mu\\[(\\d+),\\s*(\\d+)\\]$", "\\2", mu_cols))
  mu_post <- matrix(NA_real_, nrow = K, ncol = D)
  for (c in seq_along(mu_cols)) mu_post[k_idx[c], d_idx[c]] <- mean(samples[, mu_cols[c]])
  Sigma_post <- vector("list", K)
  for (k in 1:K) {
    idx <- grep(paste0("^Prec\\[[0-9]+,\\s*[0-9]+,\\s*", k, "\\]$"), colnames(samples))
    Prec_k <- matrix(colMeans(samples[, idx, drop = FALSE]), nrow = D, ncol = D, byrow = FALSE)
    S <- solve(Prec_k); S <- (S + t(S))/2
    Sigma_post[[k]] <- S
  }
  list(mu = mu_post, Sigma = Sigma_post)
}
ellipse_points <- function(mu, Sigma, r, n = 200) {
  ang <- seq(0, 2*pi, length.out = n)
  circle <- cbind(cos(ang), sin(ang))
  R <- chol(Sigma); pts <- circle %*% t(R) * r
  sweep(pts, 2, mu, FUN = "+")
}
make_aplats <- function(mu_list, Sigma_list,
                        probs = c(0.60,0.75,0.85,0.90,0.95,0.98,0.995),
                        alphas = seq(0.25,0.15,length.out=7),
                        fill_col = "grey20"){
  r_levels <- sqrt(qchisq(probs, df = 2))
  out <- vector("list", length(mu_list) * length(r_levels)); idx <- 1L
  for (k in seq_along(mu_list)) for (i in seq_along(r_levels)) {
    pts <- ellipse_points(mu_list[[k]], Sigma_list[[k]], r_levels[i])
    out[[idx]] <- data.frame(x = pts[,1], y = pts[,2],
                             comp = factor(k), level = i, alpha = alphas[i], fill = fill_col)
    idx <- idx + 1L
  }
  do.call(rbind, out)
}
mix_density_df <- function(mu_mat2d, Sigma_list2d, pi_vec, xy_grid){
  Kloc <- length(Sigma_list2d)
  dens <- rep(0, nrow(xy_grid))
  Xmat <- as.matrix(xy_grid)
  for(k in 1:Kloc){
    dens <- dens + pi_vec[k] * mvtnorm::dmvnorm(Xmat, mean = mu_mat2d[k,], sigma = Sigma_list2d[[k]], log = FALSE)
  }
  cbind(xy_grid, dens = dens)
}
get_pi_means <- function(samples, K){
  pi_cols <- grep("^pi\\[", colnames(samples), value = TRUE)
  ord <- order(as.integer(sub("^pi\\[(\\d+)\\]$", "\\1", pi_cols)))
  colMeans(samples[, pi_cols[ord], drop = FALSE])[1:K]
}
fmt_pi <- function(v) paste0("[", paste(sprintf("%.2f", v), collapse = ", "), "]")

# Emprise commune pour les plots
pad <- 1
xmin_all <- min(c(X_su_all[,1], X_wi_all[,1])) - pad
xmax_all <- max(c(X_su_all[,1], X_wi_all[,1])) + pad
ymin_all <- min(c(X_su_all[,2], X_wi_all[,2])) - pad
ymax_all <- max(c(X_su_all[,2], X_wi_all[,2])) + pad

theme0 <- theme_minimal(base_size = 12)
square_theme <- theme0 + theme(aspect.ratio = 1)
cols_season <- c("été" = "orange", "hiver" = "darkblue")

# 1) AVANT Ω (points combinés)
before_df <- rbind(
  data.frame(x = X_su_all[,1], y = X_su_all[,2], saison = "été"),
  data.frame(x = X_wi_all[,1], y = X_wi_all[,2], saison = "hiver")
)
p_before <- ggplot(before_df, aes(x = x, y = y, color = saison)) +
  geom_point(size = 0.55, alpha = 0.85) +
  scale_color_manual(values = cols_season) +
  coord_equal(xlim = c(xmin_all, xmax_all), ylim = c(ymin_all, ymax_all), expand = FALSE) +
  square_theme + labs(title = "1) AVANT Ω — Été & Hiver", x = "x", y = "y", color = "Saison")

# 2) FONDS Ω (été et hiver côte à côte)
su_r_df <- as.data.frame(r_su, xy = TRUE); names(su_r_df)[3] <- "omega"
wi_r_df <- as.data.frame(r_wi, xy = TRUE); names(wi_r_df)[3] <- "omega"
p_su_omega <- ggplot(su_r_df, aes(x = x, y = y, fill = omega)) +
  geom_raster() +
  coord_equal(xlim = c(xmin_all, xmax_all), ylim = c(ymin_all, ymax_all), expand = FALSE) +
  scale_fill_viridis_c(name = "Ω (été)") +
  square_theme + labs(title = "2) Fond Ω — ÉTÉ", x = "x", y = "y")
p_wi_omega <- ggplot(wi_r_df, aes(x = x, y = y, fill = omega)) +
  geom_raster() +
  coord_equal(xlim = c(xmin_all, xmax_all), ylim = c(ymin_all, ymax_all), expand = FALSE) +
  scale_fill_viridis_c(name = "Ω (hiver)") +
  square_theme + labs(title = "2) Fond Ω — HIVER", x = "x", y = "y")
p_on_omega <- (p_su_omega | p_wi_omega)

# 3) APRÈS Ω (points combinés)
X_su_obs <- X[,1:2,drop=FALSE]; X_wi_obs <- X[,3:4,drop=FALSE]
after_df <- rbind(
  data.frame(x = X_su_obs[,1], y = X_su_obs[,2], saison = "été"),
  data.frame(x = X_wi_obs[,1], y = X_wi_obs[,2], saison = "hiver")
)
p_after <- ggplot(after_df, aes(x = x, y = y, color = saison)) +
  geom_point(size = 0.55, alpha = 0.9) +
  scale_color_manual(values = cols_season) +
  coord_equal(xlim = c(xmin_all, xmax_all), ylim = c(ymin_all, ymax_all), expand = FALSE) +
  square_theme + labs(title = "3) APRÈS Ω — Été & Hiver", x = "x", y = "y", color = "Saison")

# ---------- Projections 4D -> 2D pour densités/ellipses ----------
post_naif <- extr_muSigma4D(samples_naif, K, D = 4)
post_corr <- extr_muSigma4D(samples_corr, K, D = 4)

mu_true_mat <- do.call(rbind, mu)
mu_su_true  <- mu_true_mat[, 1:2, drop=FALSE]
mu_wi_true  <- mu_true_mat[, 3:4, drop=FALSE]
Sig_su_true <- lapply(true_Sigma, \(S) S[1:2, 1:2, drop=FALSE])
Sig_wi_true <- lapply(true_Sigma, \(S) S[3:4, 3:4, drop=FALSE])

mu_su_naif  <- post_naif$mu[, 1:2, drop=FALSE]
mu_wi_naif  <- post_naif$mu[, 3:4, drop=FALSE]
Sig_su_naif <- lapply(post_naif$Sigma, \(S) S[1:2, 1:2, drop=FALSE])
Sig_wi_naif <- lapply(post_naif$Sigma, \(S) S[3:4, 3:4, drop=FALSE])

mu_su_corr  <- post_corr$mu[, 1:2, drop=FALSE]
mu_wi_corr  <- post_corr$mu[, 3:4, drop=FALSE]
Sig_su_corr <- lapply(post_corr$Sigma, \(S) S[1:2, 1:2, drop=FALSE])
Sig_wi_corr <- lapply(post_corr$Sigma, \(S) S[3:4, 3:4, drop=FALSE])

su_xy <- as.data.frame(r_su, xy = TRUE)[, c("x","y")]
wi_xy <- as.data.frame(r_wi, xy = TRUE)[, c("x","y")]

dens_su_true <- as.data.frame(mix_density_df(mu_su_true, Sig_su_true, pi_true, su_xy))
dens_su_naif <- as.data.frame(mix_density_df(mu_su_naif, Sig_su_naif, get_pi_means(samples_naif, K), su_xy))
dens_su_corr <- as.data.frame(mix_density_df(mu_su_corr, Sig_su_corr, get_pi_means(samples_corr, K), su_xy))

dens_wi_true <- as.data.frame(mix_density_df(mu_wi_true, Sig_wi_true, pi_true, wi_xy))
dens_wi_naif <- as.data.frame(mix_density_df(mu_wi_naif, Sig_wi_naif, get_pi_means(samples_naif, K), wi_xy))
dens_wi_corr <- as.data.frame(mix_density_df(mu_wi_corr, Sig_wi_corr, get_pi_means(samples_corr, K), wi_xy))

mk_polys <- function(mu2d_mat, Sig2d_list, col){
  polys <- make_aplats(split(mu2d_mat, row(mu2d_mat)), Sig2d_list, fill_col = col)
  polys[order(polys$level, decreasing = TRUE), ]
}
polys_su_true <- mk_polys(mu_su_true, Sig_su_true, "grey20")
polys_su_naif <- mk_polys(mu_su_naif, Sig_su_naif, "darkorange3")
polys_su_corr <- mk_polys(mu_su_corr, Sig_su_corr, "blue4")

polys_wi_true <- mk_polys(mu_wi_true, Sig_wi_true, "grey20")
polys_wi_naif <- mk_polys(mu_wi_naif, Sig_wi_naif, "darkorange3")
polys_wi_corr <- mk_polys(mu_wi_corr, Sig_wi_corr, "blue4")

df_su_true_cent <- data.frame(x = mu_su_true[,1], y = mu_su_true[,2])
df_su_naif_cent <- data.frame(x = mu_su_naif[,1], y = mu_su_naif[,2])
df_su_corr_cent <- data.frame(x = mu_su_corr[,1], y = mu_su_corr[,2])

df_wi_true_cent <- data.frame(x = mu_wi_true[,1], y = mu_wi_true[,2])
df_wi_naif_cent <- data.frame(x = mu_wi_naif[,1], y = mu_wi_naif[,2])
df_wi_corr_cent <- data.frame(x = mu_wi_corr[,1], y = mu_wi_corr[,2])

p_gmm <- function(dens_df, centers_df, title_txt, subtitle_txt, xlim, ylim){
  ggplot() +
    geom_raster(data = dens_df, aes(x = x, y = y, fill = dens), interpolate = TRUE) +
    geom_contour(data = dens_df, aes(x = x, y = y, z = dens),
                 color = "white", alpha = 0.6, bins = 12, linewidth = 0.25) +
    geom_point(data = centers_df, aes(x = x, y = y),
               color = "black", size = 2.0, shape = 3, stroke = 0.8) +
    coord_equal(xlim = xlim, ylim = ylim, expand = FALSE) +
    scale_fill_viridis_c(name = "densité",
                         guide = guide_colorbar(barheight = unit(70, "pt"))) +
    theme0 + labs(title = title_txt, subtitle = subtitle_txt, x = "x", y = "y")
}
p_ell <- function(polys, cents, col, xlim, ylim){
  ggplot() +
    geom_polygon(data = polys,
                 aes(x = x, y = y, group = interaction(comp, level), alpha = alpha),
                 fill = col) + scale_alpha_identity() +
    geom_point(data = cents, aes(x = x, y = y),
               color = col,
               shape = if (col=="darkorange3") 17 else if (col=="blue4") 18 else 4,
               size = 3.0, stroke = if (col=="grey20") 1.0 else 0.8) +
    coord_equal(xlim = xlim, ylim = ylim, expand = FALSE) +
    theme0 + labs(x = "x", y = "y")
}

xlim_su <- c(xmin_all, xmax_all); ylim_su <- c(ymin_all, ymax_all)
xlim_wi <- c(xmin_all, xmax_all); ylim_wi <- c(ymin_all, ymax_all)

pi_naif <- get_pi_means(samples_naif, K)
pi_corr <- get_pi_means(samples_corr, K)
pi_true <- pi_true

# Été — densités + ellipses
p_su_d_true <- p_gmm(dens_su_true, df_su_true_cent, "ÉTÉ — Densité GMM (vrai)",    paste0("π = ", fmt_pi(pi_true)), xlim_su, ylim_su)
p_su_d_naif <- p_gmm(dens_su_naif, df_su_naif_cent, "ÉTÉ — Densité GMM (naïf)",    paste0("π = ", fmt_pi(pi_naif)), xlim_su, ylim_su)
p_su_d_corr <- p_gmm(dens_su_corr, df_su_corr_cent, "ÉTÉ — Densité GMM (corrigé)", paste0("π = ", fmt_pi(pi_corr)), xlim_su, ylim_su)

p_su_e_true <- p_ell(polys_su_true, df_su_true_cent, "grey20",     xlim_su, ylim_su)
p_su_e_naif <- p_ell(polys_su_naif, df_su_naif_cent, "darkorange3",xlim_su, ylim_su)
p_su_e_corr <- p_ell(polys_su_corr, df_su_corr_cent, "blue4",      xlim_su, ylim_su)

# Hiver — densités + ellipses
p_wi_d_true <- p_gmm(dens_wi_true, df_wi_true_cent, "HIVER — Densité GMM (vrai)",    paste0("π = ", fmt_pi(pi_true)), xlim_wi, ylim_wi)
p_wi_d_naif <- p_gmm(dens_wi_naif, df_wi_naif_cent, "HIVER — Densité GMM (naïf)",    paste0("π = ", fmt_pi(pi_naif)), xlim_wi, ylim_wi)
p_wi_d_corr <- p_gmm(dens_wi_corr, df_wi_corr_cent, "HIVER — Densité GMM (corrigé)", paste0("π = ", fmt_pi(pi_corr)), xlim_wi, ylim_wi)

p_wi_e_true <- p_ell(polys_wi_true, df_wi_true_cent, "grey20",     xlim_wi, ylim_wi)
p_wi_e_naif <- p_ell(polys_wi_naif, df_wi_naif_cent, "darkorange3",xlim_wi, ylim_wi)
p_wi_e_corr <- p_ell(polys_wi_corr, df_wi_corr_cent, "blue4",      xlim_wi, ylim_wi)

# ---------- Assemblages ----------



row_top <- (p_before | p_after) / (p_su_omega | p_wi_omega)

lab_panel <- function(txt){
  ggplot() + theme_void() + coord_cartesian(clip = "off") +
    annotate("text", x = 1, y = 0.5, label = txt,
             angle = 90, fontface = "bold", size = 4.2, hjust = 1) +
    theme(plot.margin = margin(0, -12, 0, 0))
}
lab_density  <- lab_panel("Densité GMM")
lab_ellipses <- lab_panel("Ellipses")
labels_width <- 0.08

su_cols <- ((p_su_d_true / p_su_e_true) | (p_su_d_naif / p_su_e_naif) | (p_su_d_corr / p_su_e_corr))
wi_cols <- ((p_wi_d_true / p_wi_e_true) | (p_wi_d_naif / p_wi_e_naif) | (p_wi_d_corr / p_wi_e_corr))


(su_cols) / (wi_cols) 
