# ==================================================
# GMM Bayésien (NIMBLE) + raster Ω (vectorisé) + rejet selon Ω
# Comparaison "corrigé par Ω" vs "naïf (ignore Ω)"
# ==================================================
# 1) Simulation des deux ellipses réelles
# 2) Écriture d’Ω
# 3) Application d’Ω => points observés
# 4) GMM Bayésien avec prise en compte d’Ω (ones-trick)
# 5) GMM Bayésien naïf (ignore Ω)
# 6) Plots : ellipses estimées (corrigé vs naïf) + ellipses réelles
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
# 1) Simulation des données (2D, Σ pleines)
# --------------------------------------------------
N_tot <- 2000
D <- 2
K <- 2

true_mu <- list(c(-2, 0),
                c( 0.5, 1))
true_sigma <- list(
  matrix(c(0.8,  0.4,
           0.4,  0.6), 2, 2, byrow = TRUE),
  matrix(c(1.5, -0.5,
           -0.5,  1.2), 2, 2, byrow = TRUE)
)
true_pi <- c(0.30, 0.55)

z <- sample(1:K, N_tot, replace = TRUE, prob = true_pi)
X_all <- t(vapply(
  1:N_tot,
  function(i) MASS::mvrnorm(1, true_mu[[z[i]]], true_sigma[[z[i]]]),
  numeric(D)
))

# --------------------------------------------------
# 2) Raster Ω(x) vectorisé + définition d’Ω
# --------------------------------------------------
xmin <- min(X_all[,1]) - 1; xmax <- max(X_all[,1]) + 1
ymin <- min(X_all[,2]) - 1; ymax <- max(X_all[,2]) + 1

nx <- 50; ny <- 50
r <- rast(nrows = ny, ncols = nx, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)

xy_grid <- as.data.frame(xyFromCell(r, 1:ncell(r)))
colnames(xy_grid) <- c("x","y")

# Ω vectorisée
omega_fun <- function(x, y){
  base <- 0.3
  spot1 <- 3 * exp(-((x - (-2.5))^2 + (y - 2.5)^2) / 4)
  spot2 <- 0.8 * exp(-((x - (-2.5))^2 + (y - (-0.25))^2) / 4)
  left <- (x < 0)
  val <- pmin(1, base + spot1)
  val <- pmax(0, val - spot2)
  val[left] <- val[left] * 0.30
  val
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

# (Objets utiles pour Z si besoin)
A_cell <- prod(res(r))
grid <- as.matrix(xy_grid)
omega_grid <- as.numeric(omega_vals)

# --------------------------------------------------
# 4) Modèle GMM (corrigé par Ω) | NIMBLE + ones-trick
# --------------------------------------------------
# Densité MVN stable via Cholesky en nimbleFunction
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
    ll[i] <- log(Omega[i] * mixdens[i] + 1e-300)  # <= 0 en pratique ici
    ones[i] ~ dbern(exp(ll[i]))                   # ones-trick
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

samples_corr <- runMCMC(cmcmc_corr, niter = 4000, nburnin = 2000, thin = 4, setSeed = TRUE)

# --------------------------------------------------
# 5) Modèle GMM naïf (ignore Ω) | NIMBLE (z latent)
# --------------------------------------------------
# Modèle standard : z[i] ~ Cat(pi),  X[i,] ~ N(mu[z[i]], Prec[,,z[i]])
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
data_naif <- list(
  X = X
)
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

samples_naif <- runMCMC(cmcmc_naif, niter = 4000, nburnin = 2000, thin = 4, setSeed = TRUE)

# --------------------------------------------------
# 6) Plots : ellipses Réelles vs Corrigé(Ω) vs Naïf
# --------------------------------------------------

# ---------- outils ellipses ----------
ellipse_points <- function(mu, Sigma, r, n = 200) {
  ang <- seq(0, 2*pi, length.out = n)
  circle <- cbind(cos(ang), sin(ang))
  R <- chol(Sigma)
  pts <- circle %*% t(R) * r
  sweep(pts, 2, mu, FUN = "+")
}
make_aplats <- function(mu_list, Sigma_list,
                        probs = c(0.60, 0.75, 0.85, 0.90, 0.95, 0.98, 0.995),
                        alphas = seq(0.25, 0.15, length.out = 7),
                        fill_col = "grey20") {
  r_levels <- sqrt(qchisq(probs, df = 2))
  out <- vector("list", length(mu_list) * length(r_levels))
  idx <- 1L
  for (k in seq_along(mu_list)) {
    for (i in seq_along(r_levels)) {
      pts <- ellipse_points(mu_list[[k]], Sigma_list[[k]], r_levels[i])
      out[[idx]] <- data.frame(
        x = pts[,1], y = pts[,2],
        comp = factor(k),
        level = i,
        alpha = alphas[i],
        fill = fill_col
      )
      idx <- idx + 1L
    }
  }
  do.call(rbind, out)
}

# ---------- post-traitement : paramètres corrigés ----------
mu_cols_c <- grep("^mu\\[", colnames(samples_corr), value = TRUE)
k_idx_c <- as.integer(sub("^mu\\[(\\d+),\\s*(\\d+)\\]$", "\\1", mu_cols_c))
d_idx_c <- as.integer(sub("^mu\\[(\\d+),\\s*(\\d+)\\]$", "\\2", mu_cols_c))
K_c <- max(k_idx_c); D_c <- max(d_idx_c)
stopifnot(K_c == K, D_c == D)

mu_post_c <- matrix(NA_real_, nrow = K, ncol = D)
for (c in seq_along(mu_cols_c)) {
  mu_post_c[k_idx_c[c], d_idx_c[c]] <- mean(samples_corr[, mu_cols_c[c]])
}
sigma_est_c <- vector("list", K)
for (k in 1:K) {
  idx <- grep(paste0("^Prec\\[[0-9]+,\\s*[0-9]+,\\s*", k, "\\]$"), colnames(samples_corr))
  prec_vals <- colMeans(samples_corr[, idx, drop = FALSE])
  Prec_k <- matrix(prec_vals, nrow = D, ncol = D, byrow = FALSE)
  Sigma_k <- solve(Prec_k); Sigma_k <- (Sigma_k + t(Sigma_k))/2
  sigma_est_c[[k]] <- Sigma_k
}

# ---------- post-traitement : paramètres naïfs ----------
mu_cols_n <- grep("^mu\\[", colnames(samples_naif), value = TRUE)
k_idx_n <- as.integer(sub("^mu\\[(\\d+),\\s*(\\d+)\\]$", "\\1", mu_cols_n))
d_idx_n <- as.integer(sub("^mu\\[(\\d+),\\s*(\\d+)\\]$", "\\2", mu_cols_n))

mu_post_n <- matrix(NA_real_, nrow = K, ncol = D)
for (c in seq_along(mu_cols_n)) {
  mu_post_n[k_idx_n[c], d_idx_n[c]] <- mean(samples_naif[, mu_cols_n[c]])
}
sigma_est_n <- vector("list", K)
for (k in 1:K) {
  idx <- grep(paste0("^Prec\\[[0-9]+,\\s*[0-9]+,\\s*", k, "\\]$"), colnames(samples_naif))
  prec_vals <- colMeans(samples_naif[, idx, drop = FALSE])
  Prec_k <- matrix(prec_vals, nrow = D, ncol = D, byrow = FALSE)
  Sigma_k <- solve(Prec_k); Sigma_k <- (Sigma_k + t(Sigma_k))/2
  sigma_est_n[[k]] <- Sigma_k
}

# ---------- préparer aplats ----------
polys_true <- make_aplats(true_mu, true_sigma, fill_col = "grey20")
polys_cor  <- make_aplats(split(mu_post_c, row(mu_post_c)),
                          sigma_est_c, fill_col = "blue4")
polys_nai  <- make_aplats(split(mu_post_n, row(mu_post_n)),
                          sigma_est_n, fill_col = "darkorange3")

# dessiner de l’extérieur vers l’intérieur
ord <- function(df) df[order(df$level, decreasing = TRUE), ]
polys_true <- ord(polys_true)
polys_cor  <- ord(polys_cor)
polys_nai  <- ord(polys_nai)

# centres
df_true_cent <- data.frame(x = do.call(rbind, true_mu)[,1],
                           y = do.call(rbind, true_mu)[,2], what = "Vrai")
df_cor_cent  <- data.frame(x = mu_post_c[,1],
                           y = mu_post_c[,2], what = "Corrigé Ω")
df_nai_cent  <- data.frame(x = mu_post_n[,1],
                           y = mu_post_n[,2], what = "Naïf")

# données observées + fond Ω
r_df <- as.data.frame(r, xy = TRUE); names(r_df)[3] <- "omega"
df_obs <- data.frame(x = X[,1], y = X[,2])

p_omega <- ggplot(r_df, aes(x = x, y = y, fill = omega)) +
  geom_raster() +
  coord_equal() +
  scale_fill_viridis_c(name = "Ω", limits = c(0,1)) +
  theme_minimal(base_size = 12) +
  labs(title = "Fond Ω et points observés") +
  geom_point(data = df_obs, aes(x = x, y = y), inherit.aes = FALSE, size = 0.6)

# plot principal : ellipses
p_ellipses <- ggplot() +
  # vrais
  geom_polygon(data = polys_true,
               aes(x = x, y = y, group = interaction(comp, level), alpha = alpha),
               fill = "grey20") +
  # corrigé Ω
  geom_polygon(data = polys_cor,
               aes(x = x, y = y, group = interaction(comp, level), alpha = alpha),
               fill = "blue4") +
  # naïf
  geom_polygon(data = polys_nai,
               aes(x = x, y = y, group = interaction(comp, level), alpha = alpha),
               fill = "darkorange3") +
  scale_alpha_identity() +
  geom_point(data = df_obs, aes(x = x, y = y),
             color = "black", alpha = 0.8, size = 0.5) +
  geom_point(data = df_true_cent, aes(x = x, y = y),
             color = "red", shape = 4, size = 3.6, stroke = 1.1) +
  geom_point(data = df_cor_cent, aes(x = x, y = y),
             color = "blue4", shape = 18, size = 3.4) +
  geom_point(data = df_nai_cent, aes(x = x, y = y),
             color = "darkorange3", shape = 17, size = 3.2) +
  coord_equal() +
  theme_minimal(base_size = 13) +
  labs(title = "Ellipses : Réelles (gris), Corrigé Ω (bleu), Naïf (orange)",
       x = "x", y = "y")

# Affichage combiné
(p_omega / p_ellipses) + plot_layout(heights = c(1, 1.2))
