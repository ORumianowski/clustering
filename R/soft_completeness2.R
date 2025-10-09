# ==================================================
# GMM Bayésien (NIMBLE) + raster Ω (vectorisé) + rejet selon Ω
# Ones-trick dans le modèle, -N log Z calculé en post-traitement
# ==================================================
library(nimble)
library(MASS)
library(terra)
library(mvtnorm)
library(ggplot2)
library(patchwork)

set.seed(123)

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

df <- as.data.frame(X_all)
colnames(df) <- c("x1", "x2")
p_X_all =  ggplot(df, aes(x = x1, y = x2)) +
  geom_point(color = "blue", size = 1) +
  theme_minimal()

# --------------------------------------------------
# 2) Raster Ω(x) vectorisé + rejet Bernoulli
# --------------------------------------------------
xmin <- min(X_all[,1]) - 1; xmax <- max(X_all[,1]) + 1
ymin <- min(X_all[,2]) - 1; ymax <- max(X_all[,2]) + 1

nx <- 50; ny <- 50
r <- rast(nrows = ny, ncols = nx, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)

xy_grid <- as.data.frame(xyFromCell(r, 1:ncell(r)))
colnames(xy_grid) <- c("x","y")

# >>> omega_fun VECTORISÉE (celle d'avant) <<<
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

Omega_all <- terra::extract(r, X_all, method = "simple")[,1]
keep <- runif(N_tot) < Omega_all
X <- X_all[keep, , drop = FALSE]
N <- nrow(X)
Omega_vec <- Omega_all[keep]
cat("Nombre de points observés après Ω :", N, "\n")

df <- as.data.frame(X)
colnames(df) <- c("x1", "x2")
p_X = ggplot(df, aes(x = x1, y = x2)) +
  geom_point(color = "blue", size = 1) +
  theme_minimal()

# === A) Préparer un ggplot du raster Ω ===
# -> convertir le raster terra en data.frame pour ggplot
r_df <- as.data.frame(r, xy = TRUE)
names(r_df)[3] <- "omega"

p_omega <- ggplot(r_df, aes(x = x, y = y, fill = omega)) +
  geom_raster() +
  coord_equal() +
  scale_fill_viridis_c(name = "Ω", limits = c(0, 1)) +
  theme_minimal() +
  labs(title = "Raster Ω(x)")

# === B) Superposer les points après Ω sur le fond Ω ===
# (on utilise les points X = après rejet)
df_after <- as.data.frame(X)
colnames(df_after) <- c("x1", "x2")

p_omega_after <- p_omega +
  geom_point(data = df_after, aes(x = x1, y = x2),
             inherit.aes = FALSE, size = 0.7) +
  labs(title = "Points retenus sur fond Ω")

# (Optionnel) superposer AVANT + APRÈS avec des couleurs différentes :
df_before <- as.data.frame(X_all); colnames(df_before) <- c("x1", "x2")
p_omega_before_after <- p_omega +
  geom_point(data = df_before, aes(x = x1, y = x2),
             inherit.aes = FALSE, size = 0.6, alpha = 0.25) +
  geom_point(data = df_after,  aes(x = x1, y = x2),
             inherit.aes = FALSE, size = 0.7) +
  labs(title = "Avant (transparent) + Après (plein) sur Ω")

# === C) Composition des figures ===
# 1) Avant vs Après (comme tu l’avais : p_X_all + p_X)
fig_av_ap <- p_X_all + p_X

# 2) Plot combiné avec Ω (choisis l’un des deux ci-dessous)
fig_omega_apres <- p_omega_after
# ou, pour voir avant+après en même temps :
fig_omega_avant_apres <- p_omega_before_after

# 3) Tout ensemble en une figure : avant/après en haut, Ω+points en bas
(fig_av_ap / fig_omega_apres) +
  plot_layout(heights = c(1, 1.15))


# Constantes pour Z
## Aire d'un pixel
A_cell <- prod(res(r))
## Nombre de pixels
M <- ncell(r)
# Coordonnées des pixels
grid <- as.matrix(xy_grid)
# Valeurs des pixels (d'effort d'échantillonnage)
omega_grid <- as.numeric(omega_vals)

# --------------------------------------------------
# 3) Fonction R pour log Z(θ,Ω) (post-traitement)
# --------------------------------------------------
# A partir de omega et des paramtres des guassiennes: calcule log(Z)

compute_logZ <- function(pi, mu, Prec){
  dens_mix <- rowSums(sapply(1:K, function(k){
    mvtnorm::dmvnorm(grid, mean = mu[k,], sigma = solve(Prec[,,k]), log = FALSE) * pi[k]
  }))
  Z <- A_cell * sum(omega_grid * dens_mix)
  log(pmax(Z, 1e-300))
}

# --------------------------------------------------
# 4) Densité MVN en nimbleFunction (stable via Cholesky)
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

# --------------------------------------------------
# 5) Modèle NIMBLE avec ones-trick (sans -logZ dans le modèle)
# --------------------------------------------------
# Remarque : on met p_i = exp(ll_i) dans un Bernoulli. Pour rester dans [0,1],
# on s'appuie sur le fait que log(Omega * densité) <= 0 avec ces Σ (densité <~ 0.3).
code <- nimbleCode({
  for(i in 1:N){
    for(k in 1:K){
      dens[i,k] <- pi[k] * exp(dmvnorm_nimble(X[i,1:D], mu[k,1:D], Prec[1:D,1:D,k], TRUE))
    }
    mixdens[i] <- sum(dens[i,1:K])
    ll[i] <- log(Omega[i] * mixdens[i] + 1e-300)  # <= 0 dans notre paramétrage
    ones[i] ~ dbern(exp(ll[i]))                   # ones-trick
  }
  
  pi[1:K] ~ ddirch(alpha[1:K])
  for(k in 1:K){
    mu[k,1:D] ~ dmnorm(mu0[1:D], prec = Prec0[1:D,1:D])
    Prec[1:D,1:D,k] ~ dwish(R[1:D,1:D], df)
  }
})

# --------------------------------------------------
# 6) Constantes / données / inits
# --------------------------------------------------
Prec0 <- diag(D) * 1e-2
R <- diag(D)

constants <- list(
  N = N, D = D, K = K,
  alpha = rep(1, K),
  mu0 = c(0, 0),
  Prec0 = Prec0,
  R = R,
  df = D + 2
)
data <- list(
  X = X,
  Omega = Omega_vec,
  ones = rep(1, N)
)

Prec_init <- array(0, dim = c(D, D, K))
for(k in 1:K) Prec_init[,,k] <- diag(D)

inits <- list(
  mu  = matrix(rnorm(K*D), K, D),
  Prec = Prec_init,
  pi  = rep(1/K, K)
)

# --------------------------------------------------
# 7) Compilation & MCMC
# --------------------------------------------------
model  <- nimbleModel(code, data = data, constants = constants, inits = inits)
cmodel <- compileNimble(model)

conf <- configureMCMC(model, monitors = c("mu","pi","Prec"))
mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc, project = model)

samples <- runMCMC(cmcmc, niter = 3000, nburnin = 1500, thin = 3)
print(head(samples))


library(ggplot2)
library(ggplot2)
library(mvtnorm)

library(ggplot2)

# ---------- 1) barycentres estimés ----------
mu_cols <- grep("^mu\\[", colnames(samples), value = TRUE)
k_idx <- as.integer(sub("^mu\\[(\\d+),\\s*(\\d+)\\]$", "\\1", mu_cols))
d_idx <- as.integer(sub("^mu\\[(\\d+),\\s*(\\d+)\\]$", "\\2", mu_cols))
K <- max(k_idx); D <- max(d_idx)

mu_post <- matrix(NA_real_, nrow = K, ncol = D)
for (c in seq_along(mu_cols)) {
  mu_post[k_idx[c], d_idx[c]] <- mean(samples[, mu_cols[c]])
}

# ---------- 2) covariances estimées ----------
sigma_est <- vector("list", K)
for (k in 1:K) {
  idx <- grep(paste0("^Prec\\[[0-9]+,\\s*[0-9]+,\\s*", k, "\\]$"), colnames(samples))
  prec_vals <- colMeans(samples[, idx, drop = FALSE])
  Prec_k <- matrix(prec_vals, nrow = D, ncol = D, byrow = FALSE)
  Sigma_k <- solve(Prec_k)
  Sigma_k <- (Sigma_k + t(Sigma_k)) / 2  # symétriser
  sigma_est[[k]] <- Sigma_k
}

# ---------- 3) fonction pour aplats ----------
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
  out <- list()
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

# ---------- 4) aplats vrais + estimés ----------
polys_true <- make_aplats(true_mu, true_sigma, fill_col = "grey20")
polys_est  <- make_aplats(split(mu_post, row(mu_post)),
                          sigma_est, fill_col = "blue4")

# ordonner pour dessiner de l’extérieur vers l’intérieur
polys_true <- polys_true[order(polys_true$level, decreasing = TRUE), ]
polys_est  <- polys_est[order(polys_est$level,  decreasing = TRUE), ]

# ---------- 5) centres ----------
df_data <- data.frame(x = X[,1], y = X[,2])
df_true_cent <- data.frame(x = do.call(rbind, true_mu)[,1],
                           y = do.call(rbind, true_mu)[,2], what = "Vrai")
df_est_cent  <- data.frame(x = mu_post[,1],
                           y = mu_post[,2], what = "Estimé")

# ---------- 6) plot ----------
p <- ggplot() +
  # aplats vrais (gris)
  geom_polygon(data = polys_true,
               aes(x = x, y = y, group = interaction(comp, level), alpha = alpha),
               fill = "grey20") +
  # aplats estimés (bleu)
  geom_polygon(data = polys_est,
               aes(x = x, y = y, group = interaction(comp, level), alpha = alpha),
               fill = "blue4") +
  scale_alpha_identity() +
  # données observées
  geom_point(data = df_data, aes(x = x, y = y),
             color = "black", alpha = 0.9, size = 0.5) +
  # centres vrais et estimés
  geom_point(data = df_true_cent, aes(x = x, y = y),
             color = "red", shape = 4, size = 3.8, stroke = 1.2) +
  geom_point(data = df_est_cent, aes(x = x, y = y),
             color = "blue", shape = 18, size = 3.8) +
  coord_equal() +
  theme_minimal(base_size = 13) +
  labs(title = "Isoclines en aplats : Vrais (gris) et Estimés (bleu)",
       x = "x", y = "y")

print(p)

library(MASS)
library(ggplot2)

set.seed(123)

# --------------------------------------------------
# 1. Simulation des données
# --------------------------------------------------
N_all <- 400        # nombre de points simulés avant rejet
D <- 2; K <- 2

true_mu <- list(c(-2, 0), c(3, 3))
true_sigma <- list(diag(2)*0.5, diag(2))
true_pi <- c(0.4, 0.6)

# tirer les composantes
z <- rbinom(N_all, 1, true_pi[2]) + 1

X_all <- t(sapply(1:N_all, function(i) {
  if(z[i]==1) MASS::mvrnorm(1, true_mu[[1]], true_sigma[[1]])
  else MASS::mvrnorm(1, true_mu[[2]], true_sigma[[2]])
}))

# fonction effort Ω(x)
effort_fun <- function(x) {
  if(x[1] < 0) return(0.3) else return(1.0)
}

Omega_all <- apply(X_all, 1, effort_fun)

# tirage Bernoulli pour savoir si un point est retenu
is_sampled <- rbinom(N_all, 1, Omega_all) == 1

X <- X_all[is_sampled, , drop = FALSE]   # données observées

# --------------------------------------------------
# 2. Exemple: après MCMC, on a 'samples' (déjà calculé dans ton script)
# Ici je suppose que 'samples' existe déjà
# --------------------------------------------------

# ---------- barycentres estimés ----------
mu_cols <- grep("^mu\\[", colnames(samples), value = TRUE)
k_idx <- as.integer(sub("^mu\\[(\\d+),\\s*(\\d+)\\]$", "\\1", mu_cols))
d_idx <- as.integer(sub("^mu\\[(\\d+),\\s*(\\d+)\\]$", "\\2", mu_cols))
K <- max(k_idx); D <- max(d_idx)

mu_post <- matrix(NA_real_, nrow = K, ncol = D)
for (c in seq_along(mu_cols)) {
  mu_post[k_idx[c], d_idx[c]] <- mean(samples[, mu_cols[c]])
}

# ---------- covariances estimées ----------
sigma_est <- vector("list", K)
for (k in 1:K) {
  idx <- grep(paste0("^Prec\\[[0-9]+,\\s*[0-9]+,\\s*", k, "\\]$"),
              colnames(samples))
  prec_vals <- colMeans(samples[, idx, drop = FALSE])
  Prec_k <- matrix(prec_vals, nrow = D, ncol = D, byrow = FALSE)
  Sigma_k <- solve(Prec_k)
  Sigma_k <- (Sigma_k + t(Sigma_k)) / 2
  sigma_est[[k]] <- Sigma_k
}

# ---------- fonctions pour aplats ----------
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
  out <- list()
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

# ---------- aplats ----------
polys_true <- make_aplats(true_mu, true_sigma, fill_col = "grey20")
polys_est  <- make_aplats(split(mu_post, row(mu_post)),
                          sigma_est, fill_col = "blue4")

polys_true <- polys_true[order(polys_true$level, decreasing = TRUE), ]
polys_est  <- polys_est[order(polys_est$level,  decreasing = TRUE), ]

# ---------- centres ----------
df_true_cent <- data.frame(x = do.call(rbind, true_mu)[,1],
                           y = do.call(rbind, true_mu)[,2], what = "Vrai")
df_est_cent  <- data.frame(x = mu_post[,1],
                           y = mu_post[,2], what = "Estimé")

# ---------- données avec statut échantillonnage ----------
df_data <- data.frame(
  x = X_all[,1],
  y = X_all[,2],
  sampled = factor(is_sampled,
                   levels = c(FALSE, TRUE),
                   labels = c("non échantillonné", "échantillonné"))
)

# --------------------------------------------------
# 3. Plot final
# --------------------------------------------------
p <- ggplot() +
  # aplats vrais (gris)
  geom_polygon(data = polys_true,
               aes(x = x, y = y, group = interaction(comp, level), alpha = alpha),
               fill = "grey20") +
  # aplats estimés (bleu)
  geom_polygon(data = polys_est,
               aes(x = x, y = y, group = interaction(comp, level), alpha = alpha),
               fill = "blue4") +
  scale_alpha_identity() +
  # points échantillonnés / non
  geom_point(data = df_data, aes(x = x, y = y, color = sampled),
             alpha = 0.7, size = 0.9) +
  scale_color_manual(values = c("non échantillonné" = "grey70",
                                "échantillonné"    = "black")) +
  # centres
  geom_point(data = df_true_cent, aes(x = x, y = y),
             color = "red", shape = 4, size = 3.8, stroke = 1.2) +
  geom_point(data = df_est_cent, aes(x = x, y = y),
             color = "blue", shape = 18, size = 3.8) +
  coord_equal() +
  theme_minimal(base_size = 13) +
  labs(title = "Isoclines : Vrais (gris) et Estimés (bleu)\nPoints : échantillonnés (noir) vs non (gris)",
       x = "x", y = "y", color = "Statut")

print(p)
