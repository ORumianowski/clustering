# ======================== PLOTS (nouvelle version) ========================
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(mvtnorm); library(viridisLite) })

# ---------- Helpers ----------
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
make_aplats_sets <- function(mu_list_su, Sigma_list_su,
                             mu_list_wi, Sigma_list_wi,
                             probs = c(0.50, 0.75, 0.90)) {
  r_levels   <- sqrt(qchisq(probs, df = 2))
  cov_labels <- sprintf("%d%%", round(100*probs))
  out <- vector("list", 2 * length(mu_list_su) * length(r_levels))
  idx <- 1L
  for (k in seq_along(mu_list_su)) {
    # été
    for (i in seq_along(r_levels)) {
      pts <- ellipse_points(mu_list_su[[k]], Sigma_list_su[[k]], r_levels[i])
      out[[idx]] <- data.frame(
        x = pts[,1], y = pts[,2],
        cluster  = factor(k, levels = as.character(1:K)),
        coverage = factor(cov_labels[i], levels = c("50%","75%","90%")),
        set = factor("été", levels = c("été","hiver"))
      ); idx <- idx + 1L
    }
    # hiver
    for (i in seq_along(r_levels)) {
      pts <- ellipse_points(mu_list_wi[[k]], Sigma_list_wi[[k]], r_levels[i])
      out[[idx]] <- data.frame(
        x = pts[,1], y = pts[,2],
        cluster  = factor(k, levels = as.character(1:K)),
        coverage = factor(cov_labels[i], levels = c("50%","75%","90%")),
        set = factor("hiver", levels = c("été","hiver"))
      ); idx <- idx + 1L
    }
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
assign_clusters_MAP <- function(X4, mu_mat, Sigma_list, pi_vec){
  K <- nrow(mu_mat)
  dens <- sapply(1:K, function(k) pi_vec[k] * mvtnorm::dmvnorm(X4, mean = mu_mat[k,], sigma = Sigma_list[[k]], log = FALSE))
  max.col(dens, ties.method = "first")
}

# ---------- Emprise / thèmes ----------
pad <- 1
xmin_all <- min(c(X_su_all[,1], X_wi_all[,1])) - pad
xmax_all <- max(c(X_su_all[,1], X_wi_all[,1])) + pad
ymin_all <- min(c(X_su_all[,2], X_wi_all[,2])) - pad
ymax_all <- max(c(X_su_all[,2], X_wi_all[,2])) + pad
xlim_all <- c(xmin_all, xmax_all); ylim_all <- c(ymin_all, ymax_all)

theme0 <- theme_minimal(base_size = 12)
square_theme <- theme0 + theme(aspect.ratio = 1)

# ---------- Densités (inchangées) ÉTÉ / HIVER ----------
post_naif <- extr_muSigma4D(samples_naif, K, D = 4)
post_corr <- extr_muSigma4D(samples_corr, K, D = 4)

# Vrais paramètres projetés
mu_true_mat <- do.call(rbind, mu)
mu_su_true  <- mu_true_mat[, 1:2, drop=FALSE]
mu_wi_true  <- mu_true_mat[, 3:4, drop=FALSE]
Sig_su_true <- lapply(true_Sigma, \(S) S[1:2, 1:2, drop=FALSE])
Sig_wi_true <- lapply(true_Sigma, \(S) S[3:4, 3:4, drop=FALSE])

# Naïf projeté
mu_su_naif  <- post_naif$mu[, 1:2, drop=FALSE]
mu_wi_naif  <- post_naif$mu[, 3:4, drop=FALSE]
Sig_su_naif <- lapply(post_naif$Sigma, \(S) S[1:2, 1:2, drop=FALSE])
Sig_wi_naif <- lapply(post_naif$Sigma, \(S) S[3:4, 3:4, drop=FALSE])

# Corrigé projeté
mu_su_corr  <- post_corr$mu[, 1:2, drop=FALSE]
mu_wi_corr  <- post_corr$mu[, 3:4, drop=FALSE]
Sig_su_corr <- lapply(post_corr$Sigma, \(S) S[1:2, 1:2, drop=FALSE])
Sig_wi_corr <- lapply(post_corr$Sigma, \(S) S[3:4, 3:4, drop=FALSE])

# Grilles XY à partir des rasters
su_xy <- as.data.frame(r_su, xy = TRUE)[, c("x","y")]
wi_xy <- as.data.frame(r_wi, xy = TRUE)[, c("x","y")]

dens_su_true <- as.data.frame(mix_density_df(mu_su_true, Sig_su_true, pi_true, su_xy))
dens_su_naif <- as.data.frame(mix_density_df(mu_su_naif, Sig_su_naif, get_pi_means(samples_naif, K), su_xy))
dens_su_corr <- as.data.frame(mix_density_df(mu_su_corr, Sig_su_corr, get_pi_means(samples_corr, K), su_xy))

dens_wi_true <- as.data.frame(mix_density_df(mu_wi_true, Sig_wi_true, pi_true, wi_xy))
dens_wi_naif <- as.data.frame(mix_density_df(mu_wi_naif, Sig_wi_naif, get_pi_means(samples_naif, K), wi_xy))
dens_wi_corr <- as.data.frame(mix_density_df(mu_wi_corr, Sig_wi_corr, get_pi_means(samples_corr, K), wi_xy))

p_gmm <- function(dens_df, centers_df, title_txt, subtitle_txt){
  ggplot() +
    geom_raster(data = dens_df, aes(x = x, y = y, fill = dens), interpolate = TRUE) +
    geom_contour(data = dens_df, aes(x = x, y = y, z = dens),
                 color = "white", alpha = 0.6, bins = 12, linewidth = 0.25) +
    geom_point(data = centers_df, aes(x = x, y = y),
               color = "black", size = 2.0, shape = 3, stroke = 0.8) +
    coord_equal(xlim = xlim_all, ylim = ylim_all, expand = FALSE) +
    scale_fill_viridis_c(name = "densité") +
    theme0 + labs(title = title_txt, subtitle = subtitle_txt, x = "x", y = "y")
}

df_su_true_cent <- data.frame(x = mu_su_true[,1], y = mu_su_true[,2])
df_su_naif_cent <- data.frame(x = mu_su_naif[,1], y = mu_su_naif[,2])
df_su_corr_cent <- data.frame(x = mu_su_corr[,1], y = mu_su_corr[,2])

df_wi_true_cent <- data.frame(x = mu_wi_true[,1], y = mu_wi_true[,2])
df_wi_naif_cent <- data.frame(x = mu_wi_naif[,1], y = mu_wi_naif[,2])
df_wi_corr_cent <- data.frame(x = mu_wi_corr[,1], y = mu_wi_corr[,2])

pi_naif <- get_pi_means(samples_naif, K)
pi_corr <- get_pi_means(samples_corr, K)

p_su_d_true <- p_gmm(dens_su_true, df_su_true_cent, "ÉTÉ — Densité GMM (vrai)",    paste0("π = [", paste(sprintf("%.2f", pi_true), collapse = ", "), "]"))
p_su_d_naif <- p_gmm(dens_su_naif, df_su_naif_cent, "ÉTÉ — Densité GMM (naïf)",    paste0("π = [", paste(sprintf("%.2f", pi_naif), collapse = ", "), "]"))
p_su_d_corr <- p_gmm(dens_su_corr, df_su_corr_cent, "ÉTÉ — Densité GMM (corrigé)", paste0("π = [", paste(sprintf("%.2f", pi_corr), collapse = ", "), "]"))

p_wi_d_true <- p_gmm(dens_wi_true, df_wi_true_cent, "HIVER — Densité GMM (vrai)",    paste0("π = [", paste(sprintf("%.2f", pi_true), collapse = ", "), "]"))
p_wi_d_naif <- p_gmm(dens_wi_naif, df_wi_naif_cent, "HIVER — Densité GMM (naïf)",    paste0("π = [", paste(sprintf("%.2f", pi_naif), collapse = ", "), "]"))
p_wi_d_corr <- p_gmm(dens_wi_corr, df_wi_corr_cent, "HIVER — Densité GMM (corrigé)", paste0("π = [", paste(sprintf("%.2f", pi_corr), collapse = ", "), "]"))

# ---------- Ellipses combinées (été + hiver) + segments individuels ----------
# TRUE labels (connus avant filtrage) pour les points conservés
z_true_obs <- factor(z[keep], levels = as.character(1:K))

# MAP labels (naïf & corrigé) sur les points observés (en 4D)
assign_naif <- assign_clusters_MAP(X, post_naif$mu, post_naif$Sigma, pi_naif)
assign_corr <- assign_clusters_MAP(X, post_corr$mu, post_corr$Sigma, pi_corr)
z_naif_obs  <- factor(assign_naif, levels = as.character(1:K))
z_corr_obs  <- factor(assign_corr, levels = as.character(1:K))

# DataFrames segments été->hiver
pairs_true <- data.frame(
  lonA = X[,1], latA = X[,2], lonB = X[,3], latB = X[,4],
  cluster = z_true_obs
)
pairs_naif <- transform(pairs_true, cluster = z_naif_obs)
pairs_corr <- transform(pairs_true, cluster = z_corr_obs)

# Ellipses combinées par modèle (été + hiver sur le même plan)
mk_polys_sets <- function(mu2d_su, Sig2d_su, mu2d_wi, Sig2d_wi){
  make_aplats_sets(
    mu_list_su = split(mu2d_su, row(mu2d_su)),
    Sigma_list_su = Sig2d_su,
    mu_list_wi = split(mu2d_wi, row(mu2d_wi)),
    Sigma_list_wi = Sig2d_wi,
    probs = c(0.50, 0.75, 0.90)
  )
}
polys_true <- mk_polys_sets(mu_su_true, Sig_su_true, mu_wi_true, Sig_wi_true)
polys_naif <- mk_polys_sets(mu_su_naif, Sig_su_naif, mu_wi_naif, Sig_wi_naif)
polys_corr <- mk_polys_sets(mu_su_corr, Sig_su_corr, mu_wi_corr, Sig_wi_corr)

# Palette & alphas
cluster_palette <- setNames(rocket(K, begin = 0.10, end = 0.75, direction = -1),
                            as.character(1:K))
coverage_alphas <- c("50%" = 0.50, "75%" = 0.20, "90%" = 0.10)

# Builder combiné (segments + ellipses été/hiver)
build_pairs_plot <- function(pairs_df, polys_df, title){
  ggplot() +
    # Segments été -> hiver
    geom_segment(data = pairs_df,
                 aes(x = lonA, y = latA, xend = lonB, yend = latB, color = cluster),
                 alpha = 0.13, linewidth = 0.5) +
    # Ellipses été & hiver superposées (sets)
    geom_polygon(data = polys_df,
                 aes(x = x, y = y, group = interaction(cluster, coverage, set),
                     fill = cluster, alpha = coverage),
                 color = NA) +
    scale_fill_manual(values = cluster_palette, name = "Clusters") +
    scale_color_manual(values = cluster_palette, name = "Clusters") +
    scale_alpha_manual(values = coverage_alphas, name = "Probability contours") +
    guides(color = guide_legend(order = 1),
           fill  = guide_legend(order = 2),
           alpha = guide_legend(order = 3)) +
    coord_equal(xlim = xlim_all, ylim = ylim_all, expand = FALSE) +
    theme0 + labs(title = title, x = "x", y = "y")
}

p_pairs_true <- build_pairs_plot(pairs_true, polys_true, "Ellipses combinées — VRAI (été + hiver)")
p_pairs_naif <- build_pairs_plot(pairs_naif, polys_naif, "Ellipses combinées — NAÏF (été + hiver)")
p_pairs_corr <- build_pairs_plot(pairs_corr, polys_corr, "Ellipses combinées — CORRIGÉ (été + hiver)")

# ---------- Assemblage final ----------
# ---------- 1) Densités : deux rangées (été puis hiver) ----------
row_su <- (p_su_d_true | p_su_d_naif | p_su_d_corr)
row_wi <- (p_wi_d_true | p_wi_d_naif | p_wi_d_corr)

p_top <- (row_su / row_wi) +
  plot_annotation(title = "Densités séparées — Été (haut) / Hiver (bas)") &
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Afficher la figure des densités
print(p_top)

# ---------- 2) Ellipses combinées + segments : 1 seule légende ----------
row_pairs <- (p_pairs_true | p_pairs_naif | p_pairs_corr) +
  plot_layout(guides = "collect") &             # <-- collecte les légendes
  theme(legend.position = "bottom")              # <-- place la légende unique en bas

# Afficher la figure ellipses+segments (légende unique)
print(row_pairs)
