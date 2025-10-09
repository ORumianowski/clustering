# ==================================================
# GMM Bayésien (NIMBLE) avec correction Ω :
# log p_o(x|θ,Ω) = log Ω(x) + log p(x|θ) - log Z(θ,Ω)
# Évalué via ones-trick: 1 ~ Bernoulli(exp(ℓ_i - C_ub))
# (C_ub grand, constant et indépendant de θ)
# ==================================================

rm(list = ls())
set.seed(123)

suppressPackageStartupMessages({
  library(nimble)
  library(MASS)
  library(terra)
  library(mvtnorm)
  library(ggplot2)
  library(patchwork)
})

# ----------------------- utilitaires -----------------------
gauss2d <- function(x, y, mu, Sigma) {
  invS <- solve(Sigma)
  dx <- x - mu[1]; dy <- y - mu[2]
  Q <- invS[1,1]*dx*dx + (invS[1,2]+invS[2,1])*dx*dy + invS[2,2]*dy*dy
  exp(-0.5 * Q)
}

ellipse_points <- function(mu, Sigma, r, n = 200) {
  ang <- seq(0, 2*pi, length.out = n)
  circle <- cbind(cos(ang), sin(ang))
  R <- chol(Sigma); pts <- circle %*% t(R) * r
  sweep(pts, 2, mu, FUN = "+")
}
make_aplats <- function(mu_list, Sigma_list,
                        probs = c(0.60, 0.75, 0.85, 0.90, 0.95, 0.98, 0.995),
                        alphas = seq(0.25, 0.15, length.out = 7),
                        fill_col = "grey20") {
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
extr_muSigma <- function(samples, K, D){
  mu_cols <- grep("^mu\\[", colnames(samples), value = TRUE)
  k_idx <- as.integer(sub("^mu\\[(\\d+),\\s*(\\d+)\\]$", "\\1", mu_cols))
  d_idx <- as.integer(sub("^mu\\[(\\d+),\\s*(\\d+)\\]$", "\\2", mu_cols))
  mu_post <- matrix(NA_real_, nrow = K, ncol = D)
  for (c in seq_along(mu_cols)) mu_post[k_idx[c], d_idx[c]] <- mean(samples[, mu_cols[c]])
  sigma_est <- vector("list", K)
  for (k in 1:K) {
    idx <- grep(paste0("^Prec\\[[0-9]+,\\s*[0-9]+,\\s*", k, "\\]$"), colnames(samples))
    Prec_k <- matrix(colMeans(samples[, idx, drop = FALSE]), nrow = D, ncol = D, byrow = FALSE)
    Sigma_k <- solve(Prec_k); Sigma_k <- (Sigma_k + t(Sigma_k))/2
    sigma_est[[k]] <- Sigma_k
  }
  list(mu = mu_post, Sigma = sigma_est)
}

# ----------------------- 1) Simulation -----------------------
N_tot <- 5000; D <- 2; K <- 3

true_mu <- list(
  c(-3.2,  0.0),   # C1 : vertical gauche
  c( 0.2, -0.6),   # C2 : diagonale (-45°)
  c( 3.2,  0.0)    # C3 : vertical droite
)

# C1 et C3 : verticales ; C2 : diagonale (-45°) via Sigma explicite
true_sigma <- list(
  # C1 : verticale gauche (x un peu plus large, y nettement plus large)
  matrix(c(0.80, 0.00,
           0.00, 4.00), 2, 2, byrow = TRUE),
  
  # C2 : diagonale (-45°) plus dispersée
  # Formule explicite pour angle -45° : 
  # diag = (a+b)/2, off = -(a-b)/2, avec a=10.0 (major), b=1.2 (minor)
  matrix(c(5.60, -4.40,
           -4.40,  5.60), 2, 2, byrow = TRUE),
  
  # C3 : verticale droite (x un peu plus large, y large)
  matrix(c(0.90, 0.00,
           0.00, 3.60), 2, 2, byrow = TRUE)
)


true_pi <- c(0.20, 0.45, 0.35)  # C2 majoritaire

z <- sample(1:K, N_tot, TRUE, true_pi)
X_all <- t(vapply(1:N_tot, function(i) MASS::mvrnorm(1,true_mu[[z[i]]],true_sigma[[z[i]]]), numeric(D)))

# ----------------------- 2) Ω(x) indépendante -----------------------
xmin <- min(X_all[,1])-1; xmax <- max(X_all[,1])+1
ymin <- min(X_all[,2])-1; ymax <- max(X_all[,2])+1
nx <- 60; ny <- 60
r <- rast(nrows=ny,ncols=nx,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax)
xy_grid <- as.data.frame(xyFromCell(r,1:ncell(r))); colnames(xy_grid) <- c("x","y")

omega_fun <- function(x, y){
  base <- 0.04
  # Bosses hautes (centres/Σ/amp FIXES — indépendants des vraies gaussiennes)
  bumpL <- 0.90 * gauss2d(x, y, c(-3.8, 2.6), diag(c(0.60, 0.60)))
  bumpR <- 0.90 * gauss2d(x, y, c( 3.8, 2.6), diag(c(0.65, 0.65)))
  # Grand trou bas horizontal (indépendant)
  bumpC <- 1.00 * gauss2d(x, y, c( 0.0, -3.2), diag(c(6.0, 0.9)))
  pmax(0, pmin(base + bumpL  + bumpR + bumpC , 0.70))
}

omega_vals <- with(xy_grid, omega_fun(x,y))
r[] <- omega_vals

# ----------------------- 3) Rejet Ω -----------------------
Omega_all <- terra::extract(r, X_all, method="simple")[,1]
keep <- runif(N_tot) < Omega_all
X <- X_all[keep,,drop=FALSE]; N <- nrow(X)
Omega_vec <- Omega_all[keep]
cat("N observés après Ω :", N, "\n")

# ----------------------- 4) Grille pour Z(θ,Ω) -----------------------
A_cell     <- prod(res(r))
M          <- ncell(r)
grid       <- as.matrix(xy_grid[,c("x","y")])
omega_grid <- as.numeric(omega_vals)

# ----------------------- 5) Modèle corrigé Ω -----------------------
# ====== 5) Modèle corrigé Ω : ones-trick stable et borné ======

# MVN (Cholesky)
dmvnorm_nimble <- nimbleFunction(
  run = function(x = double(1), mean = double(1), Prec = double(2),
                 log = logical(0, default = FALSE)) {
    returnType(double(0))
    Dloc <- length(x); xm <- x - mean
    qf <- inprod(xm, Prec %*% xm)
    U <- chol(Prec); ldet <- 2 * sum(log(diag(U)))
    logdens <- 0.5 * ldet - 0.5 * Dloc * log(2*pi) - 0.5 * qf
    if (log) return(logdens) else return(exp(logdens))
  }
)

# log Z(θ,Ω) par quadrature (grille fixe)
logZ_calc <- nimbleFunction(
  run = function(mu   = double(2),   # K x D
                 Prec = double(3),   # D x D x K
                 pi   = double(1),   # K
                 grid = double(2),   # M x D
                 omega= double(1),   # M
                 A    = double(0),
                 K    = integer(0),
                 D    = integer(0),
                 M    = integer(0)) {
    returnType(double(0))
    sumZ <- 0.0
    for(m in 1:M){
      mix <- 0.0
      for(k in 1:K){
        mix <- mix + pi[k] * dmvnorm_nimble(grid[m,1:D], mu[k,1:D],
                                            Prec[1:D,1:D,k], FALSE)
      }
      sumZ <- sumZ + omega[m] * mix
    }
    Z <- A * sumZ
    if (Z < 1e-300) Z <- 1e-300
    return(log(Z))
  })

# --- constantes numériques pour le ones-trick ---
C_ub   <- 200       # assez grand, constant et indépendant de θ
p_min  <- 1e-300    # bornes de sécurité [p_min, 1 - p_min]
p_max1 <- 1 - 1e-12

code_corrige <- nimbleCode({
  # terme commun à tous : log Z(θ,Ω)
  logZ <- logZ_calc(mu[1:K,1:D], Prec[1:D,1:D,1:K], pi[1:K],
                    grid[1:M,1:D], omega[1:M], A, K, D, M)
  
  for (i in 1:N) {
    # mélange p(x_i | θ)
    for (k in 1:K) {
      dens[i,k] <- pi[k] * exp(dmvnorm_nimble(X[i,1:D], mu[k,1:D],
                                              Prec[1:D,1:D,k], TRUE))
    }
    mixdens[i] <- sum(dens[i,1:K])
    
    # log po(x_i | θ,Ω) = log Ω(x_i) + log mixdens - logZ
    ll[i] <- log(Omega[i] * mixdens[i] + 1e-300) - logZ
    
    # ones-trick: p_i = exp(ll - C_ub) CLIPPE dans (0,1)
    p_raw[i]      <- exp(ll[i] - C_ub)
    p_clip_hi[i]  <- min(p_raw[i], p_max1)
    p[i]          <- max(p_clip_hi[i], p_min)
    ones[i] ~ dbern(p[i])
  }
  
  # Priors
  pi[1:K] ~ ddirch(alpha[1:K])
  for (k in 1:K) {
    mu[k,1:D] ~ dmnorm(mu0[1:D], prec = Prec0[1:D,1:D])
    Prec[1:D,1:D,k] ~ dwish(R[1:D,1:D], df)
  }
})

# --- constants / data pour ce modèle (si pas déjà faits) ---
Prec0 <- diag(D) * 1e-2
R     <- diag(D)

constants_corr <- list(
  N = N, D = D, K = K,
  M = M, A = A_cell,
  alpha = rep(1, K),
  mu0 = c(0, 0),
  Prec0 = Prec0,
  R = R, df = D + 2,
  # bornes du ones-trick (scalaires constants)
  p_min = 1e-300, p_max1 = 1 - 1e-12, C_ub = 200
)
data_corr <- list(
  X = X,
  Omega = Omega_vec,
  ones = rep(1, N),
  grid = grid,
  omega = omega_grid
)

# inits
Prec_init <- array(0, dim = c(D, D, K)); for(k in 1:K) Prec_init[,,k] <- diag(D)
inits_corr <- list(mu = matrix(rnorm(K*D), K, D), Prec = Prec_init, pi = rep(1/K, K))

# build & run
model_corr  <- nimbleModel(code_corrige, data = data_corr, constants = constants_corr, inits = inits_corr, check = FALSE)
cmodel_corr <- compileNimble(model_corr)
conf_corr   <- configureMCMC(model_corr, monitors = c("mu","pi","Prec"))
mcmc_corr   <- buildMCMC(conf_corr)
cmcmc_corr  <- compileNimble(mcmc_corr, project = model_corr)
samples_corr <- runMCMC(cmcmc_corr, niter = 2000, nburnin = 1000, thin = 8, setSeed = TRUE)

# ----------------------- 6) (Optionnel) Naïf pour comparaison -----------------------
code_naif <- nimbleCode({
  for (i in 1:N){
    z[i] ~ dcat(pi[1:K])
    X[i,1:D] ~ dmnorm(mu[z[i],1:D], prec = Prec[1:D,1:D,z[i]])
  }
  pi[1:K] ~ ddirch(alpha[1:K])
  for (k in 1:K){
    mu[k,1:D] ~ dmnorm(mu0[1:D], prec = Prec0[1:D,1:D])
    Prec[1:D,1:D,k] ~ dwish(R[1:D,1:D], df)
  }
})
constants_naif <- list(N=N,D=D,K=K,alpha=rep(1,K),mu0=c(0,0),Prec0=Prec0,R=R,df=D+2)
data_naif <- list(X=X)
inits_naif <- list(z=sample(1:K,N,TRUE), mu=matrix(rnorm(K*D),K,D), Prec=Prec_init, pi=rep(1/K,K))
model_naif  <- nimbleModel(code_naif, data=data_naif, constants=constants_naif, inits=inits_naif, check=FALSE)
cmodel_naif <- compileNimble(model_naif)
conf_naif   <- configureMCMC(model_naif, monitors=c("mu","pi","Prec"))
mcmc_naif   <- buildMCMC(conf_naif)
cmcmc_naif  <- compileNimble(mcmc_naif, project=model_naif)
samples_naif <- runMCMC(cmcmc_naif, niter=4000, nburnin=2000, thin=8, setSeed=TRUE)

# ----------------------- 7) 6 plots 

# ----- outils ellipses -----
ellipse_points <- function(mu, Sigma, r, n = 200) {
  ang <- seq(0, 2*pi, length.out = n)
  circle <- cbind(cos(ang), sin(ang))
  R <- chol(Sigma); pts <- circle %*% t(R) * r
  sweep(pts, 2, mu, FUN = "+")
}
make_aplats <- function(mu_list, Sigma_list,
                        probs = c(0.60, 0.75, 0.85, 0.90, 0.95, 0.98, 0.995),
                        alphas = seq(0.25, 0.15, length.out = 7),
                        fill_col = "grey20") {
  r_levels <- sqrt(qchisq(probs, df = 2))
  out <- vector("list", length(mu_list) * length(r_levels)); idx <- 1L
  for (k in seq_along(mu_list)) for (i in seq_along(r_levels)) {
    pts <- ellipse_points(mu_list[[k]], Sigma_list[[k]], r_levels[i])
    out[[idx]] <- data.frame(
      x = pts[,1], y = pts[,2],
      comp = factor(k), level = i, alpha = alphas[i], fill = fill_col
    ); idx <- idx + 1L
  }
  do.call(rbind, out)
}

# ----- post-traitement : paramètres corrigés -----
extr_muSigma <- function(samples, K, D){
  mu_cols <- grep("^mu\\[", colnames(samples), value = TRUE)
  k_idx <- as.integer(sub("^mu\\[(\\d+),\\s*(\\d+)\\]$", "\\1", mu_cols))
  d_idx <- as.integer(sub("^mu\\[(\\d+),\\s*(\\d+)\\]$", "\\2", mu_cols))
  mu_post <- matrix(NA_real_, nrow = K, ncol = D)
  for (c in seq_along(mu_cols)) mu_post[k_idx[c], d_idx[c]] <- mean(samples[, mu_cols[c]])
  sigma_est <- vector("list", K)
  for (k in 1:K) {
    idx <- grep(paste0("^Prec\\[[0-9]+,\\s*[0-9]+,\\s*", k, "\\]$"), colnames(samples))
    Prec_k <- matrix(colMeans(samples[, idx, drop = FALSE]), nrow = D, ncol = D, byrow = FALSE)
    Sigma_k <- solve(Prec_k); Sigma_k <- (Sigma_k + t(Sigma_k))/2
    sigma_est[[k]] <- Sigma_k
  }
  list(mu = mu_post, Sigma = sigma_est)
}
post_c <- extr_muSigma(samples_corr, K, D)
post_n <- extr_muSigma(samples_naif, K, D)

# ====== Bloc de remplacement (avec π affichés dans les 3 graphes du bas) ======
# ======================= FUSION EN UN SEUL BLOC =======================
# Construit :
#  - p1,p2,p3 : (1) Avant Ω | (2) Fond Ω | (3) Après Ω
#  - p4,p5,p6 : Densité GMM (Réel | Naïf | Corrigé) avec π en sous-titre
#  - p7,p8,p9 : Ellipses (Réel | Naïf | Corrigé)
#  - Mise en page finale : ligne du haut avec titre commun,
#    puis deux lignes fusionnées par colonne avec labels à gauche
#    et colorbars de densité à droite de CHAQUE panneau de densité
# ======================================================================
get_pi_means <- function(samples, K){
  pi_cols <- grep("^pi\\[", colnames(samples), value = TRUE)
  ord <- order(as.integer(sub("^pi\\[(\\d+)\\]$", "\\1", pi_cols)))
  colMeans(samples[, pi_cols[ord], drop = FALSE])[1:K]
}
fmt_pi <- function(v) paste0("[", paste(sprintf("%.2f", v), collapse = ", "), "]")

# --- π (réel/naïf/corrigé)
K <- length(true_mu)
pi_true <- true_pi
pi_naif <- get_pi_means(samples_naif, K)
pi_corr <- get_pi_means(samples_corr, K)

# --- Ellipses (aplats)
polys_true <- make_aplats(true_mu, true_sigma, fill_col = "grey20")
polys_cor  <- make_aplats(split(post_c$mu, row(post_c$mu)), post_c$Sigma, fill_col = "blue4")
polys_nai  <- make_aplats(split(post_n$mu, row(post_n$mu)), post_n$Sigma, fill_col = "darkorange3")
ord_ap <- function(df) df[order(df$level, decreasing = TRUE), ]
polys_true <- ord_ap(polys_true); polys_cor <- ord_ap(polys_cor); polys_nai <- ord_ap(polys_nai)

# --- Centres
df_true_cent <- data.frame(x = do.call(rbind, true_mu)[,1], y = do.call(rbind, true_mu)[,2])
df_cor_cent  <- data.frame(x = post_c$mu[,1], y = post_c$mu[,2])
df_nai_cent  <- data.frame(x = post_n$mu[,1], y = post_n$mu[,2])

# --- Données (ligne 1)
df_before <- data.frame(x = X_all[,1], y = X_all[,2])
df_after  <- data.frame(x = X[,1],    y = X[,2])
r_df <- as.data.frame(r, xy = TRUE); names(r_df)[3] <- "omega"

p_before <- ggplot(df_before, aes(x = x, y = y)) +
  geom_point(size = 0.55, alpha = 0.8, color = "steelblue4") +
  coord_equal(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
  theme_minimal(base_size = 12) + labs(title = "1) Points AVANT Ω", x = "x", y = "y")

p_on_omega <- ggplot(r_df, aes(x = x, y = y, fill = omega)) +
  geom_raster() +
  coord_equal(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
  scale_fill_viridis_c(name = "Ω", limits = c(0,1)) +
  theme_minimal(base_size = 12) + labs(title = "2) Fond Ω", x = "x", y = "y")

p_after <- ggplot(df_after, aes(x = x, y = y)) +
  geom_point(size = 0.55, alpha = 0.9, color = "black") +
  coord_equal(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
  theme_minimal(base_size = 12) + labs(title = "3) Points APRÈS Ω", x = "x", y = "y")

# --- Densité GMM (UNE densité par panneau)
mix_density_df <- function(mu_mat, Sigma_list, pi_vec, xy_grid){
  Kloc <- length(Sigma_list)
  dens <- rep(0, nrow(xy_grid))
  Xmat <- as.matrix(xy_grid)
  for(k in 1:Kloc){
    dens <- dens + pi_vec[k] * mvtnorm::dmvnorm(Xmat, mean = mu_mat[k,], sigma = Sigma_list[[k]], log = FALSE)
  }
  cbind(xy_grid, dens = dens)
}
p_gmm <- function(dens_df, centers_df, title_txt, subtitle_txt){
  ggplot() +
    geom_raster(data = dens_df, aes(x = x, y = y, fill = dens), interpolate = TRUE) +
    geom_contour(data = dens_df, aes(x = x, y = y, z = dens),
                 color = "white", alpha = 0.6, bins = 12, linewidth = 0.25,
                 inherit.aes = FALSE) +
    geom_point(data = centers_df, aes(x = x, y = y),
               inherit.aes = FALSE, color = "black", size = 2.0, shape = 3, stroke = 0.8) +
    coord_equal(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
    scale_fill_viridis_c(name = "densité",
                         guide = guide_colorbar(barheight = unit(70, "pt"))) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right",
          legend.box.margin = margin(l = 6, r = 2, t = 2, b = 2)) +
    labs(title = title_txt, subtitle = subtitle_txt, x = "x", y = "y")
}

# --- Paramètres pour densités
mu_true_mat <- do.call(rbind, true_mu)
mu_cor_mat  <- post_c$mu
mu_nai_mat  <- post_n$mu
Sig_true    <- true_sigma
Sig_cor     <- post_c$Sigma
Sig_nai     <- post_n$Sigma

dens_true <- as.data.frame(mix_density_df(mu_true_mat, Sig_true, pi_true, xy_grid))
dens_naif <- as.data.frame(mix_density_df(mu_nai_mat, Sig_nai, pi_naif, xy_grid))
dens_cor  <- as.data.frame(mix_density_df(mu_cor_mat,  Sig_cor,  pi_corr, xy_grid))

df_true_cent_d <- data.frame(x = mu_true_mat[,1], y = mu_true_mat[,2])
df_nai_cent_d  <- data.frame(x = mu_nai_mat[,1],  y = mu_nai_mat[,2])
df_cor_cent_d  <- data.frame(x = mu_cor_mat[,1],  y = mu_cor_mat[,2])

# --- Densités GMM (4–6) : π en sous-titre
p4 <- p_gmm(dens_true, df_true_cent_d, "Réel — Densité GMM",           paste0("π = ", fmt_pi(pi_true)))
p5 <- p_gmm(dens_naif, df_nai_cent_d, "Estimé naïf — Densité GMM",     paste0("π = ", fmt_pi(pi_naif)))
p6 <- p_gmm(dens_cor,  df_cor_cent_d, "Estimé corrigé — Densité GMM",  paste0("π = ", fmt_pi(pi_corr)))

# --- Ellipses (7–9)
p7 <- ggplot() +
  geom_polygon(data = polys_true,
               aes(x = x, y = y, group = interaction(comp, level), alpha = alpha),
               fill = "grey20") + scale_alpha_identity() +
  geom_point(data = df_true_cent, aes(x = x, y = y),
             color = "red", shape = 4, size = 3.0, stroke = 1.0) +
  coord_equal(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
  theme_minimal(base_size = 12) + labs(x = "x", y = "y")

p8 <- ggplot() +
  geom_polygon(data = polys_nai,
               aes(x = x, y = y, group = interaction(comp, level), alpha = alpha),
               fill = "darkorange3") + scale_alpha_identity() +
  geom_point(data = df_nai_cent, aes(x = x, y = y),
             color = "darkorange3", shape = 17, size = 3.0) +
  coord_equal(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
  theme_minimal(base_size = 12) + labs(x = "x", y = "y")

p9 <- ggplot() +
  geom_polygon(data = polys_cor,
               aes(x = x, y = y, group = interaction(comp, level), alpha = alpha),
               fill = "blue4") + scale_alpha_identity() +
  geom_point(data = df_cor_cent, aes(x = x, y = y),
             color = "blue4", shape = 18, size = 3.2) +
  coord_equal(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
  theme_minimal(base_size = 12) + labs(x = "x", y = "y")

# --- Étiquettes de lignes (collées à droite)
lab_panel <- function(txt){
  ggplot() + theme_void() + coord_cartesian(clip = "off") +
    annotate("text", x = 1, y = 0.5, label = txt,
             angle = 90, fontface = "bold", size = 4.2, hjust = 1) +
    theme(plot.margin = margin(0, -12, 0, 0))
}
lab_density  <- lab_panel("Densité GMM")
lab_ellipses <- lab_panel("Ellipses")

# --- Colonnes fusionnées (densité / ellipses)
col_real <- (p4 / p7)
col_naif <- (p5 / p8)
col_corr <- (p6 / p9)

# --- Ligne du haut (titre commun)
row_top <- (p_before | p_on_omega | p_after) +
  plot_annotation(title = "Simulation des points") &
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# --- Colonne de gauche (labels), largeur étroite
labels_width <- 0.07
col_widths   <- c(labels_width, 1, 1, 1)
left_labels  <- (lab_density / lab_ellipses) + plot_layout(heights = c(1, 1))

# --- Bloc du bas (labels + 3 colonnes) + ajustement des espacements
estim <- (left_labels | col_real | col_naif | col_corr) +
  plot_layout(widths = col_widths, guides = "keep") &
  theme(panel.spacing.x = unit(4, "pt"))

# --- Assemblage final
(row_top / estim) +
  plot_layout(heights = c(1, 2.2), guides = "keep")
