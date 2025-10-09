# ==================================================
# GMM Bayésien (NIMBLE) + Ω indépendante (ellipses fixes)
# 3 composantes simulées :
#   - C1 : verticale gauche
#   - C2 : diagonale (-45°), covariance explicite (PAS de rot_cov)
#   - C3 : verticale droite
# 6 plots : avant Ω | après Ω | points sur Ω | réelles | corrigé | naïf
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

# --------------------------------------------------
# 1) Simulation : 3 composantes (traits bleus)
# --------------------------------------------------
N_tot <- 4500
D <- 2
K <- 3

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

z <- sample(1:K, N_tot, replace = TRUE, prob = true_pi)
X_all <- t(vapply(1:N_tot, function(i)
  MASS::mvrnorm(1, true_mu[[z[i]]], true_sigma[[z[i]]]),
  numeric(D)))

# --------------------------------------------------
# 2) Ω(x) indépendante des paramètres simulés
#    (deux bosses circulaires en HAUT + grand trou bas horizontal)
# --------------------------------------------------
xmin <- min(X_all[,1]) - 1; xmax <- max(X_all[,1]) + 1
ymin <- min(X_all[,2]) - 1; ymax <- max(X_all[,2]) + 1

nx <- 120; ny <- 120
r <- rast(nrows = ny, ncols = nx, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
xy_grid <- as.data.frame(xyFromCell(r, 1:ncell(r))); colnames(xy_grid) <- c("x","y")

gauss2d <- function(x, y, mu, Sigma) {
  invS <- solve(Sigma)
  dx <- x - mu[1]; dy <- y - mu[2]
  Q <- invS[1,1]*dx*dx + (invS[1,2]+invS[2,1])*dx*dy + invS[2,2]*dy*dy
  exp(-0.5 * Q)
}

omega_fun <- function(x, y){
  base <- 0.04
  # Bosses hautes (centres/Σ/amp FIXES — indépendants des vraies gaussiennes)
  bumpL <- 0.90 * gauss2d(x, y, c(-3.8, 2.6), diag(c(0.60, 0.60)))
  bumpR <- 0.90 * gauss2d(x, y, c( 3.8, 2.6), diag(c(0.65, 0.65)))
  # Grand trou bas horizontal (indépendant)
  bumpC <- 1.00 * gauss2d(x, y, c( 0.0, -3.2), diag(c(6.0, 0.9)))
  pmax(0, pmin(base + bumpL  + bumpR + bumpC , 0.70))
}

omega_vals <- with(xy_grid, omega_fun(x, y))
r[] <- omega_vals

# --------------------------------------------------
# 3) Application d’Ω aux points simulés -> données observées
# --------------------------------------------------
Omega_all <- terra::extract(r, X_all, method = "simple")[,1]
keep <- runif(N_tot) < Omega_all
X <- X_all[keep, , drop = FALSE]
N <- nrow(X)
Omega_vec <- Omega_all[keep]
cat("Nombre de points observés après Ω :", N, "\n")

# --------------------------------------------------
# 4) GMM CORRIGÉ par Ω (ones-trick) | NIMBLE
# --------------------------------------------------
dmvnorm_nimble <- nimbleFunction(
  run = function(x = double(1), mean = double(1), Prec = double(2),
                 log = logical(0, default = FALSE)) {
    returnType(double(0))
    Dloc <- length(x)
    xm <- x - mean
    qf <- inprod(xm, Prec %*% xm)
    U <- chol(Prec)
    ldet <- 2 * sum(log(diag(U)))
    logdens <- 0.5 * ldet - 0.5 * Dloc * log(2*pi) - 0.5 * qf
    if (log) return(logdens) else return(exp(logdens))
  }
)

code_corrige <- nimbleCode({
  for(i in 1:N){
    for(k in 1:K){
      dens[i,k] <- pi[k] * exp(dmvnorm_nimble(X[i,1:D], mu[k,1:D], Prec[1:D,1:D,k], TRUE))
    }
    mixdens[i] <- sum(dens[i,1:K])
    ll[i] <- log(Omega[i] * mixdens[i] + 1e-300)    # <= 0 car Ω<=0.7 et densités modérées
    ones[i] ~ dbern(exp(ll[i]))                     # ones-trick
  }
  pi[1:K] ~ ddirch(alpha[1:K])
  for(k in 1:K){
    mu[k,1:D] ~ dmnorm(mu0[1:D], prec = Prec0[1:D,1:D])
    Prec[1:D,1:D,k] ~ dwish(R[1:D,1:D], df)
  }
})

Prec0 <- diag(D) * 1e-2
R <- diag(D)
constants_corr <- list(
  N = N, D = D, K = K,
  alpha = rep(1, K),
  mu0 = c(0, 0),
  Prec0 = Prec0,
  R = R,
  df = D + 2
)
data_corr <- list(
  X = X,
  Omega = Omega_vec,
  ones = rep(1, N)
)

Prec_init <- array(0, dim = c(D, D, K))
for(k in 1:K) Prec_init[,,k] <- diag(D)

inits_corr <- list(
  mu  = matrix(rnorm(K*D), K, D),
  Prec = Prec_init,
  pi  = rep(1/K, K)
)

model_corr  <- nimbleModel(code_corrige, data = data_corr, constants = constants_corr, inits = inits_corr, check = FALSE)
cmodel_corr <- compileNimble(model_corr)
conf_corr <- configureMCMC(model_corr, monitors = c("mu","pi","Prec"), enableWAIC = FALSE)
mcmc_corr <- buildMCMC(conf_corr)
cmcmc_corr <- compileNimble(mcmc_corr, project = model_corr)

samples_corr <- runMCMC(cmcmc_corr, niter = 6000, nburnin = 3000, thin = 6, setSeed = TRUE)

# --------------------------------------------------
# 5) GMM NAÏF (ignore Ω) | NIMBLE avec z latent
# --------------------------------------------------
code_naif <- nimbleCode({
  for (i in 1:N){
    z[i] ~ dcat(pi[1:K])
    X[i,1:D] ~ dmnorm(mu[z[i], 1:D], prec = Prec[1:D,1:D,z[i]])
  }
  pi[1:K] ~ ddirch(alpha[1:K])
  for (k in 1:K){
    mu[k,1:D] ~ dmnorm(mu0[1:D], prec = Prec0[1:D,1:D])
    Prec[1:D,1:D,k] ~ dwish(R[1:D,1:D], df)
  }
})

constants_naif <- list(
  N = N, D = D, K = K,
  alpha = rep(1, K),
  mu0 = c(0, 0),
  Prec0 = Prec0,
  R = R,
  df = D + 2
)
data_naif <- list(X = X)
inits_naif <- list(
  z   = sample(1:K, N, replace = TRUE),
  mu  = matrix(rnorm(K*D), K, D),
  Prec = Prec_init,
  pi  = rep(1/K, K)
)

model_naif  <- nimbleModel(code_naif, data = data_naif, constants = constants_naif, inits = inits_naif, check = FALSE)
cmodel_naif <- compileNimble(model_naif)
conf_naif <- configureMCMC(model_naif, monitors = c("mu","pi","Prec"), enableWAIC = FALSE)
mcmc_naif <- buildMCMC(conf_naif)
cmcmc_naif <- compileNimble(mcmc_naif, project = model_naif)

samples_naif <- runMCMC(cmcmc_naif, niter = 6000, nburnin = 3000, thin = 6, setSeed = TRUE)

# --------------------------------------------------
# 6) 6 PLOTS
# --------------------------------------------------

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

# ----- aplats -----
polys_true <- make_aplats(true_mu, true_sigma, fill_col = "grey20")
polys_cor  <- make_aplats(split(post_c$mu, row(post_c$mu)), post_c$Sigma, fill_col = "blue4")
polys_nai  <- make_aplats(split(post_n$mu, row(post_n$mu)), post_n$Sigma, fill_col = "darkorange3")

ord <- function(df) df[order(df$level, decreasing = TRUE), ]
polys_true <- ord(polys_true); polys_cor <- ord(polys_cor); polys_nai <- ord(polys_nai)

# ----- centres -----
df_true_cent <- data.frame(x = do.call(rbind, true_mu)[,1], y = do.call(rbind, true_mu)[,2])
df_cor_cent  <- data.frame(x = post_c$mu[,1], y = post_c$mu[,2])
df_nai_cent  <- data.frame(x = post_n$mu[,1], y = post_n$mu[,2])

# ----- données -----
df_before <- data.frame(x = X_all[,1], y = X_all[,2])
df_after  <- data.frame(x = X[,1],    y = X[,2])

# ----- Plots a/b/c -----
p_before <- ggplot(df_before, aes(x = x, y = y)) +
  geom_point(size = 0.55, alpha = 0.8, color = "steelblue4") +
  coord_equal(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
  theme_minimal(base_size = 12) + labs(title = "a) Points AVANT Ω", x = "x", y = "y")

r_df <- as.data.frame(r, xy = TRUE); names(r_df)[3] <- "omega"

p_after <- ggplot(df_after, aes(x = x, y = y)) +
  geom_point(size = 0.55, alpha = 0.9, color = "black") +
  coord_equal(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
  theme_minimal(base_size = 12) + labs(title = "b) Points APRÈS Ω", x = "x", y = "y")

# ----- Plots d/e/f -----
p_true <- ggplot() +
  geom_polygon(data = polys_true,
               aes(x = x, y = y, group = interaction(comp, level), alpha = alpha),
               fill = "grey20") + scale_alpha_identity() +
  geom_point(data = df_true_cent, aes(x = x, y = y),
             color = "grey20", shape = 4, size = 3.4, stroke = 1.0) +
  coord_equal(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
  theme_minimal(base_size = 12) + labs(title = "d) Ellipses RÉELLES", x = "x", y = "y")

p_corr <- ggplot() +
  geom_polygon(data = polys_cor,
               aes(x = x, y = y, group = interaction(comp, level), alpha = alpha),
               fill = "blue4") + scale_alpha_identity() +
  geom_point(data = df_cor_cent, aes(x = x, y = y),
             color = "blue4", shape = 18, size = 3.2) +
  coord_equal(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
  theme_minimal(base_size = 12) + labs(title = "e) Ellipses ESTIMÉES (corrigé Ω)", x = "x", y = "y")

p_naif <- ggplot() +
  geom_polygon(data = polys_nai,
               aes(x = x, y = y, group = interaction(comp, level), alpha = alpha),
               fill = "darkorange3") + scale_alpha_identity() +
  geom_point(data = df_nai_cent, aes(x = x, y = y),
             color = "darkorange3", shape = 17, size = 3.0) +
  coord_equal(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
  theme_minimal(base_size = 12) + labs(title = "f) Ellipses ESTIMÉES (naïf)", x = "x", y = "y")

# --- (b) Fond Ω (SANS points)
r_df <- as.data.frame(r, xy = TRUE); names(r_df)[3] <- "omega"
p_on_omega <- ggplot(r_df, aes(x = x, y = y, fill = omega)) +
  geom_raster() +
  coord_equal(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
  scale_fill_viridis_c(name = "Ω", limits = c(0,1)) +
  theme_minimal(base_size = 12) +
  labs(title = "b) Fond Ω", x = "x", y = "y")

# --- (c) Points APRÈS Ω (mettre à jour le titre si besoin)
p_after <- p_after + labs(title = "c) Points APRÈS Ω")

# --- Grille 2x3 avec les nouvelles positions :
#     Haut : (1) AVANT Ω | (2) FOND Ω | (3) APRÈS Ω
#     Bas  : (4) RÉELLES | (5) NAÏF   | (6) CORRIGÉ  [échange 5 et 6]
(p_before | p_on_omega | p_after) /
  (p_true   | p_naif     | p_corr) +
  plot_layout(guides = "collect", heights = c(1, 1.05))

