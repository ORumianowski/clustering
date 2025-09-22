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
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)

set.seed(123)

# ==========================================================
# 1) Données + standardisation SÛRE
# ==========================================================
nb_clusters <- 3
load("dataset/soft_connectivity.data.RData")

X <- cbind(
  soft_connectivity.data$longitude_breed,
  soft_connectivity.data$latitude_breed,
  soft_connectivity.data$longitude_wint,
  soft_connectivity.data$latitude_wint
) |> as.matrix()
X <- na.omit(X)
colnames(X) <- c("lon_breed","lat_breed","lon_wint","lat_wint")

safe_scale <- function(X) {
  mu <- colMeans(X)
  sdv <- apply(X, 2, sd)
  sdv[is.na(sdv) | sdv < 1e-8] <- 1.0
  Z <- sweep(X, 2, mu, "-")
  Z <- sweep(Z, 2, sdv, "/")
  attr(Z, "scaled:center") <- mu
  attr(Z, "scaled:scale")  <- sdv
  Z
}
Xz <- safe_scale(X)
X_center <- attr(Xz, "scaled:center")
X_scale  <- attr(Xz, "scaled:scale")

N <- nrow(Xz); d <- ncol(Xz); K <- nb_clusters

# ==========================================================
# 2) Modèle NIMBLE : Normal–Wishart (full covariance, très shrinké)
# ==========================================================
code_full_shrunk <- nimbleCode({
  pi[1:K] ~ ddirch(alpha[1:K])
  
  for(k in 1:K){
    # Prior sur Omega_k
    Omega[k,1:d,1:d] ~ dwish(R[1:d,1:d], df)
    
    # Precision scaled: prec_mu[k,,] = kappa0 * Omega[k,,]
    for(i in 1:d) {
      for(j in 1:d) {
        prec_mu[k,i,j] <- kappa0 * Omega[k,i,j]
      }
    }
    
    mu[k,1:d] ~ dmnorm(m0[1:d], prec = prec_mu[k,1:d,1:d])
  }
  
  for(i in 1:N){
    cl[i] ~ dcat(pi[1:K])
    X[i,1:d] ~ dmnorm(mu[cl[i],1:d], Omega[cl[i],1:d,1:d])
  }
})


# ==========================================================
# 3) Constantes + données + inits ROBUSTES
# ==========================================================
# Prior TRÈS FORT vers l'identité (diagonale) mais full-cov autorisé
R_mat   <- diag(d)         # centre du Wishart
df_wish <- d + 100         # <<< shrink puissant vers R
kappa0  <- 2.0             # <<< shrink de mu vers m0 (ici 0)

constants <- list(
  K=K, N=N, d=d,
  alpha=rep(1, K),
  m0 = rep(0, d),
  kappa0 = kappa0,
  R = R_mat,
  df = df_wish
)

data_list <- list(X = Xz)

# Inits stables
set.seed(123)
km <- kmeans(Xz, centers = K, nstart = 20)
init.mu <- km$centers  # pas ou peu de bruit au départ

# Omega init = identité (PD, logProb fini sous Wishart)
init.Omega <- array(0, dim=c(K,d,d))
for (k in 1:K) init.Omega[k,,] <- diag(d)

# cl init = kmeans
init.cl <- km$cluster
for (k in 1:K) if (!any(init.cl==k)) init.cl[sample.int(N,1)] <- k

inits_full <- list(
  cl = init.cl,
  mu = init.mu,
  Omega = init.Omega,
  pi = rep(1/K, K)
)

# ==========================================================
# 4) Helper MCMC
# ==========================================================
run_nimble <- function(code, constants, data, inits, monitors,
                       niter = 8000, nburnin = 3000, thin = 5,
                       progressBar = TRUE) {
  model  <- nimbleModel(code, constants = constants,
                        data = data, inits = inits, calculate = TRUE)
  cmodel <- compileNimble(model)
  conf   <- configureMCMC(model, monitors = monitors)
  # Conjugacy OK : dmnorm(mu|Omega) + dwish(Omega) + dmnorm(X|mu,Omega)
  mcmc   <- buildMCMC(conf)
  cmcmc  <- compileNimble(mcmc, project = model)
  samp   <- runMCMC(cmcmc, niter = niter, nburnin = nburnin,
                    thin = thin, progressBar = progressBar)
  list(model = model, samples = as.mcmc(samp))
}

# ==========================================================
# 5) Lancer le modèle
# ==========================================================
cat("Running full-covariance Gaussian mixture with strong shrinkage...\n")
res <- run_nimble(
  code_full_shrunk,
  constants,
  data_list,
  inits_full,
  monitors = c("pi","mu","Omega"),
  niter = 8000, nburnin = 3000, thin = 5
)

# ==========================================================
# 6) Échantillons -> data.frame
# ==========================================================
sdf <- as.data.frame(as.matrix(res$samples))
head(sdf, 3)

# ==========================================================
# 7) Posterior means: mu(z) et Sigma = E[Omega^{-1}] (moyenne des inverses)
# ==========================================================
# mu
mu_cols <- grep("^mu\\[", colnames(sdf), value = TRUE)
k_idx <- as.integer(sub("^mu\\[(\\d+),\\s*(\\d+)\\]$", "\\1", mu_cols))
d_idx <- as.integer(sub("^mu\\[(\\d+),\\s*(\\d+)\\]$", "\\2", mu_cols))
K <- max(k_idx); D <- max(d_idx)

mu_post_z <- matrix(NA_real_, nrow = K, ncol = D)
for (c in seq_along(mu_cols)) {
  mu_post_z[k_idx[c], d_idx[c]] <- mean(sdf[[mu_cols[c]]])
}
colnames(mu_post_z) <- colnames(X)

# Sigma(z) moyenne des inverses d'Omega
omega_name <- function(k,i,j) paste0("Omega[", k, ", ", i, ", ", j, "]")
n_samp <- nrow(sdf)
Sigma_acc <- replicate(K, matrix(0, D, D), simplify = FALSE)

for (s in 1:n_samp) {
  for (k in 1:K) {
    Om <- matrix(NA_real_, D, D)
    for (i in 1:D) for (j in 1:D) Om[i,j] <- sdf[[ omega_name(k,i,j) ]][s]
    Sig <- tryCatch(solve(Om), error=function(e) MASS::ginv(Om))
    Sigma_acc[[k]] <- Sigma_acc[[k]] + Sig
  }
}
Sigma_post_z <- lapply(Sigma_acc, function(M) M / n_samp)

# ---- Reprojection à l'échelle d'origine
unscale_mu <- function(mu_row, center, scale) center + mu_row * scale
unscale_Sigma <- function(Sigma, scale) {
  S <- diag(scale, nrow = length(scale), ncol = length(scale))
  S %*% Sigma %*% S
}
mu_post <- t(apply(mu_post_z, 1, unscale_mu, center = X_center, scale = X_scale))
colnames(mu_post) <- colnames(X)
Sigma_post <- lapply(Sigma_post_z, unscale_Sigma, scale = X_scale)

# ==========================================================
# 8) Ellipses et cartes (été = dims 1-2 ; hiver = dims 3-4)
# ==========================================================
ellipse_points <- function(mu, Sigma, r, n = 200) {
  ang <- seq(0, 2*pi, length.out = n)
  circle <- cbind(cos(ang), sin(ang))
  R <- chol(Sigma)
  pts <- circle %*% t(R) * r
  sweep(pts, 2, mu, FUN = "+")
}

make_aplats <- function(mu_list, Sigma_list,
                        probs = c(0.60, 0.75, 0.85, 0.90, 0.95),
                        alphas = seq(0.25, 0.15, length.out = 5),
                        fill_cols) {
  r_levels <- sqrt(qchisq(probs, df = 2))
  out <- list(); idx <- 1L
  for (k in seq_along(mu_list)) {
    for (i in seq_along(r_levels)) {
      pts <- ellipse_points(mu_list[[k]], Sigma_list[[k]], r_levels[i])
      out[[idx]] <- data.frame(
        x = pts[,1], y = pts[,2],
        comp = factor(k),
        level = i,
        alpha = alphas[i],
        fill = fill_cols[k]
      )
      idx <- idx + 1L
    }
  }
  dplyr::bind_rows(out) |> dplyr::arrange(desc(level))
}

cluster_cols <- c("red3", "forestgreen", "royalblue")

# Été
mu_summer    <- mu_post[, 1:2, drop = FALSE]
Sigma_summer <- lapply(Sigma_post, function(S) S[1:2, 1:2])
polys_summer <- make_aplats(split(mu_summer, row(mu_summer)),
                            Sigma_summer, fill_cols = cluster_cols)
df_cent_summer <- data.frame(x = mu_summer[,1], y = mu_summer[,2], comp = factor(1:K))

# Hiver
mu_winter    <- mu_post[, 3:4, drop = FALSE]
Sigma_winter <- lapply(Sigma_post, function(S) S[3:4, 3:4])
polys_winter <- make_aplats(split(mu_winter, row(mu_winter)),
                            Sigma_winter, fill_cols = cluster_cols)
df_cent_winter <- data.frame(x = mu_winter[,1], y = mu_winter[,2], comp = factor(1:K))

world <- ne_countries(scale = "medium", returnclass = "sf")
base_map <- geom_sf(data = world, fill = "grey92", color = "grey70", size = 0.3)

p_summer <- ggplot() +
  base_map +
  geom_polygon(data = polys_summer,
               aes(x = x, y = y, group = interaction(comp, level),
                   alpha = alpha, fill = fill), color = NA) +
  scale_alpha_identity() +
  scale_fill_identity() +
  geom_point(data = df_cent_summer, aes(x = x, y = y, color = comp),
             shape = 18, size = 3) +
  scale_color_manual(values = cluster_cols, name = "component",
                     labels = paste("Cluster", 1:K)) +
  coord_sf(xlim = c(-15, 35), ylim = c(35, 70), expand = FALSE) +
  theme_minimal(base_size = 13) +
  labs(title = "Barycentres & ellipses (full-cov, prior fort) - Été",
       x = "Longitude", y = "Latitude")

p_winter <- ggplot() +
  base_map +
  geom_polygon(data = polys_winter,
               aes(x = x, y = y, group = interaction(comp, level),
                   alpha = alpha, fill = fill), color = NA) +
  scale_alpha_identity() +
  scale_fill_identity() +
  geom_point(data = df_cent_winter, aes(x = x, y = y, color = comp),
             shape = 18, size = 3) +
  scale_color_manual(values = cluster_cols, name = "component",
                     labels = paste("Cluster", 1:K)) +
  coord_sf(xlim = c(-15, 35), ylim = c(35, 70), expand = FALSE) +
  theme_minimal(base_size = 13) +
  labs(title = "Barycentres & ellipses (full-cov, prior fort) - Hiver",
       x = "Longitude", y = "Latitude")

cowplot::plot_grid(p_summer, p_winter, ncol = 2)
