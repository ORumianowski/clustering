
# ============================
# PLOTS — 1 chaîne par modèle
# (corrigé Ω : samples_corr4 ; naïf : samples_naif4)
# ============================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
  library(stringr); library(mvtnorm); library(viridisLite); library(patchwork)
})

# --- 0) Alias des échantillons (obligatoires) ---
if (!exists("samples_corr4")) stop("samples_corr4 introuvable.")
if (!exists("samples_naif4")) stop("samples_naif4 introuvable.")
samples_corr <- as.data.frame(samples_corr4)
samples_naif <- as.data.frame(samples_naif4)

# --- 1) Helpers autonomes ---
extr_muSigma4D <- function(samples, K, D = 4){
  mu_cols <- grep("^mu\\[\\s*\\d+\\s*,\\s*\\d+\\s*\\]$", colnames(samples), value = TRUE)
  mu_post <- matrix(NA_real_, K, D)
  for (nm in mu_cols) {
    m <- regexec("^mu\\[\\s*(\\d+)\\s*,\\s*(\\d+)\\s*\\]$", nm)
    mm <- regmatches(nm, m)[[1]]; k <- as.integer(mm[2]); d <- as.integer(mm[3])
    mu_post[k, d] <- mean(samples[[nm]])
  }
  Sig <- vector("list", K)
  for (k in 1:K) {
    pat <- paste0("^Prec\\[\\s*(\\d+)\\s*,\\s*(\\d+)\\s*,\\s*", k, "\\s*\\]$")
    prec_cols <- grep(pat, colnames(samples), value = TRUE)
    if (length(prec_cols) != D*D) stop(sprintf("Colonnes Prec manquantes pour k=%d", k))
    Prec_k <- matrix(NA_real_, D, D)
    for (nm in prec_cols) {
      m <- regexec("^Prec\\[\\s*(\\d+)\\s*,\\s*(\\d+)\\s*,", nm)
      mm <- regmatches(nm, m)[[1]]; i <- as.integer(mm[2]); j <- as.integer(mm[3])
      Prec_k[i, j] <- mean(samples[[nm]])
    }
    Prec_k <- (Prec_k + t(Prec_k))/2
    ev <- eigen(Prec_k, symmetric = TRUE, only.values = TRUE)$values
    if (any(!is.finite(ev)) || any(ev <= 1e-10)) Prec_k <- Prec_k + diag(D)*1e-6
    S <- tryCatch(solve(Prec_k), error = function(e) NA)
    if (any(is.na(S))) { Prec_k <- Prec_k + diag(D)*1e-4; S <- solve(Prec_k) }
    Sig[[k]] <- (S + t(S))/2
  }
  list(mu = mu_post, Sigma = Sig)
}
get_pi_means <- function(samples, K){
  pi_cols <- grep("^pi\\[\\s*\\d+\\s*\\]$", colnames(samples), value = TRUE)
  ord <- order(as.integer(sub("^pi\\[\\s*(\\d+)\\s*\\]$", "\\1", pi_cols)))
  colMeans(samples[, pi_cols[ord], drop = FALSE])[1:K]
}
assign_MAP <- function(X4, mu_mat, Sigma_list, pi_vec){
  dens <- sapply(1:nrow(mu_mat), function(k)
    pi_vec[k] * mvtnorm::dmvnorm(X4, mean = mu_mat[k,], sigma = Sigma_list[[k]], log = FALSE))
  max.col(dens, ties.method = "first")
}
ellipse_points <- function(mu, Sigma, r, n = 200){
  ang <- seq(0, 2*pi, length.out = n); circle <- cbind(cos(ang), sin(ang))
  R <- chol(Sigma); pts <- circle %*% t(R) * r; sweep(pts, 2, mu, "+")
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
  out$level <- factor(out$level, levels=cov_labels)
  out
}
build_pairs_plot <- function(pairs_df, polys_df, title, subtitle, xlim_all, ylim_all,
                             cluster_palette, coverage_alphas){
  ggplot() +
    geom_segment(data = pairs_df,
                 aes(x=lonA, y=latA, xend=lonB, yend=latB, color=cluster),
                 alpha=.06, linewidth=.5) +
    geom_polygon(data = polys_df,
                 aes(x=x, y=y, group=interaction(cluster, level, set),
                     fill=cluster, alpha=level), color=NA) +
    scale_fill_manual(values=cluster_palette, name="Clusters") +
    scale_color_manual(values=cluster_palette, name="Clusters") +
    scale_alpha_manual(values=coverage_alphas, name="Contours prob.") +
    coord_equal(xlim=xlim_all, ylim=ylim_all, expand=FALSE) +
    theme_minimal(base_size=12) +
    labs(title=title, subtitle=subtitle, x="x", y="y")
}
fmt_pi <- function(v) paste0("π = [", paste(sprintf("%.2f", v), collapse=", "), "]")

# --- 2) Postérieurs & π (par modèle) ---
if (!exists("K")) stop("K introuvable (nombre de composantes).")
post_c <- extr_muSigma4D(samples_corr, K, D=4)
post_n <- extr_muSigma4D(samples_naif, K, D=4)
pi_c   <- get_pi_means(samples_corr, K)
pi_n   <- get_pi_means(samples_naif, K)

mu_su_c <- post_c$mu[,1:2,drop=FALSE]; mu_wi_c <- post_c$mu[,3:4,drop=FALSE]
mu_su_n <- post_n$mu[,1:2,drop=FALSE]; mu_wi_n <- post_n$mu[,3:4,drop=FALSE]
Sig_su_c <- lapply(post_c$Sigma, function(S) S[1:2,1:2,drop=FALSE])
Sig_wi_c <- lapply(post_c$Sigma, function(S) S[3:4,3:4,drop=FALSE])
Sig_su_n <- lapply(post_n$Sigma, function(S) S[1:2,1:2,drop=FALSE])
Sig_wi_n <- lapply(post_n$Sigma, function(S) S[3:4,3:4,drop=FALSE])

# --- 3) Emprise pour les panneaux 2D et segments été→hiver ---
if (!exists("X_su_all") || !exists("X_wi_all")) stop("X_su_all / X_wi_all introuvables.")
pad <- 1
xmin_all <- min(c(X_su_all[,1], X_wi_all[,1])) - pad
xmax_all <- max(c(X_su_all[,1], X_wi_all[,1])) + pad
ymin_all <- min(c(X_su_all[,2], X_wi_all[,2])) - pad
ymax_all <- max(c(X_su_all[,2], X_wi_all[,2])) + pad
xlim_all <- c(xmin_all, xmax_all); ylim_all <- c(ymin_all, ymax_all)

# --- 4) Segments été→hiver (individus observés au moins une fois) ---
if (!exists("X4_all") || !exists("sel_ete") || !exists("sel_hiver"))
  stop("X4_all / sel_ete / sel_hiver introuvables.")
X_obs_any <- X4_all[(sel_ete | sel_hiver), , drop=FALSE]
z_corr <- factor(assign_MAP(X_obs_any, post_c$mu, post_c$Sigma, pi_c))
z_naif <- factor(assign_MAP(X_obs_any, post_n$mu, post_n$Sigma, pi_n))

mk_pairs <- function(X4, zfac) data.frame(
  lonA=X4[,1], latA=X4[,2], lonB=X4[,3], latB=X4[,4], cluster=zfac
)
set.seed(1)
use_idx <- seq_len(nrow(X_obs_any))
if (length(use_idx) > 1000) use_idx <- sample(use_idx, 1000)
pairs_corr <- mk_pairs(X_obs_any, z_corr)[use_idx,]
pairs_naif <- mk_pairs(X_obs_any, z_naif)[use_idx,]

# --- 5) Ellipses (été+hiver) pour chaque modèle ---
polys_corr <- make_aplats_sets(mu_su_c, Sig_su_c, mu_wi_c, Sig_wi_c)
polys_naif <- make_aplats_sets(mu_su_n, Sig_su_n, mu_wi_n, Sig_wi_n)

# --- 6) Option “Vérité” si dispo (sinon ignoré) ---
has_truth <- FALSE
if (exists("mu_true") && exists("Sigma_true")) {
  mu_true_mat <- do.call(rbind, mu_true)
  Sigma_true_list <- Sigma_true
  has_truth <- TRUE
} else if (exists("mu_su_true") && exists("mu_wi_true") &&
           exists("Sigma_su_true") && exists("Sigma_wi_true")) {
  mu_true_mat <- do.call(rbind, lapply(1:K, function(k) c(mu_su_true[[k]], mu_wi_true)))
  Sigma_true_list <- lapply(1:K, function(k) {
    S <- matrix(0, 4, 4)
    S[1:2,1:2] <- Sigma_su_true[[k]]
    S[3:4,3:4] <- Sigma_wi_true
    S
  })
  has_truth <- TRUE
}
if (has_truth) {
  mu_su_true <- mu_true_mat[,1:2,drop=FALSE]
  mu_wi_true <- mu_true_mat[,3:4,drop=FALSE]
  Sig_su_true <- lapply(Sigma_true_list, function(S) S[1:2,1:2,drop=FALSE])
  Sig_wi_true <- lapply(Sigma_true_list, function(S) S[3:4,3:4,drop=FALSE])
  polys_true <- make_aplats_sets(mu_su_true, Sig_su_true, mu_wi_true, Sig_wi_true)
  # segments "vrais" si z existe, sinon segments neutres
  if (exists("z")) {
    z_true_obs <- factor(z[(sel_ete | sel_hiver)])
    pairs_true <- mk_pairs(X_obs_any, z_true_obs)[use_idx,]
  } else {
    pairs_true <- mk_pairs(X_obs_any, factor(1))[use_idx,]
  }
}

# --- 7) Stylisation & palettes ---
cluster_palette <- setNames(rocket(K, begin=0.10, end=0.75, direction=-1), as.character(1:K))
coverage_alphas <- c("50%"=.50, "75%"=.20, "90%"=.10)

# --- 8) Panneaux Ellipses+Segments (Naïf | Corrigé | [Vrai si possible]) ---
p_pairs_naif <- build_pairs_plot(
  pairs_naif, polys_naif,
  "Ellipses combinées — NAÏF (été + hiver)", fmt_pi(pi_n),
  xlim_all, ylim_all, cluster_palette, coverage_alphas
)
p_pairs_corr <- build_pairs_plot(
  pairs_corr, polys_corr,
  "Ellipses combinées — CORRIGÉ (été + hiver)", fmt_pi(pi_c),
  xlim_all, ylim_all, cluster_palette, coverage_alphas
)

if (has_truth) {
  p_pairs_true <- build_pairs_plot(
    pairs_true, polys_true,
    "Ellipses combinées — VRAI (été + hiver)",
    if (exists("pi_true")) fmt_pi(pi_true) else "π (simulation)",
    xlim_all, ylim_all, cluster_palette, coverage_alphas
  )
  row_pairs <- (p_pairs_true | p_pairs_naif | p_pairs_corr) +
    plot_layout(guides="collect") & theme(legend.position="bottom")
} else {
  row_pairs <- (p_pairs_naif | p_pairs_corr) +
    plot_layout(guides="collect") & theme(legend.position="bottom")
}
print(row_pairs)

# --- 9) (option) Traceplots & densités μ/π pour comparaison modèle ---
#       Facilement désactivables si tu veux seulement les ellipses/segments.
lab_dims <- c("x_su","y_su","x_wi","y_wi")
to_long <- function(S, tag){ S$.iter <- seq_len(nrow(S)); S$.model <- tag; S }
S_corr <- to_long(samples_corr, "corrigé Ω")
S_naif <- to_long(samples_naif, "naïf")
mu_cols <- grep("^mu\\[\\s*\\d+\\s*,\\s*\\d+\\s*\\]$", colnames(S_corr), value=TRUE)
pi_cols <- grep("^pi\\[\\s*\\d+\\s*\\]$", colnames(S_corr), value=TRUE)
mk_mu_map <- function(cols){
  do.call(rbind, lapply(cols, function(nm){
    m <- stringr::str_match(nm, "^mu\\[\\s*(\\d+)\\s*,\\s*(\\d+)\\s*\\]$")
    k <- as.integer(m[1,2]); d <- as.integer(m[1,3])
    data.frame(param=nm, k=k, d=d, label=sprintf("mu[k=%d, %s]", k, lab_dims[d]))
  }))
}
mu_map <- mk_mu_map(mu_cols)
long_mu <- bind_rows(
  S_corr |> select(all_of(mu_cols), .iter, .model) |> pivot_longer(all_of(mu_cols), names_to="param", values_to="value"),
  S_naif |> select(all_of(mu_cols), .iter, .model) |> pivot_longer(all_of(mu_cols), names_to="param", values_to="value")
) |> left_join(mu_map, by="param")
long_pi <- bind_rows(
  S_corr |> select(all_of(pi_cols), .iter, .model) |> pivot_longer(all_of(pi_cols), names_to="param", values_to="value"),
  S_naif |> select(all_of(pi_cols), .iter, .model) |> pivot_longer(all_of(pi_cols), names_to="param", values_to="value")
) |> mutate(k = as.integer(sub("^pi\\[(\\d+)\\]$", "\\1", param)),
            label = sprintf("pi[k=%d]", k))

p_trace_mu <- ggplot(long_mu, aes(.iter, value, color=.model)) +
  geom_line(alpha=.85, linewidth=.5) +
  facet_wrap(~ label, scales="free_y", ncol=3) +
  labs(title="Traceplots — μ (par modèle)", x="itération", y="valeur", color="modèle") +
  theme_minimal(base_size=12)

p_dens_mu <- ggplot(long_mu, aes(value, color=.model, fill=.model)) +
  geom_density(alpha=.25, linewidth=.7) +
  facet_wrap(~ label, scales="free", ncol=3) +
  labs(title="Densités postérieures — μ", x="valeur", y="densité", color="modèle", fill="modèle") +
  theme_minimal(base_size=12)

p_trace_pi <- ggplot(long_pi, aes(.iter, value, color=.model)) +
  geom_line(alpha=.85, linewidth=.6) +
  facet_wrap(~ label, scales="free_y") +
  labs(title="Traceplots — π", x="itération", y="valeur", color="modèle") +
  theme_minimal(base_size=12)

p_dens_pi <- ggplot(long_pi, aes(value, color=.model, fill=.model)) +
  geom_density(alpha=.25, linewidth=.8) +
  facet_wrap(~ label, scales="free") +
  labs(title="Densités postérieures — π", x="valeur", y="densité", color="modèle", fill="modèle") +
  theme_minimal(base_size=12)

print(p_trace_mu); print(p_dens_mu); print(p_trace_pi); print(p_dens_pi)



# ==================================================
# === FIGURE : points simulés + Ω + points observés
#      (bloC adapté au script courant)
# ==================================================

# bornes carrées
pad <- 1
xmin_all <- min(c(X_su_all[,1], X_wi_all[,1])) - pad
xmax_all <- max(c(X_su_all[,1], X_wi_all[,1])) + pad
ymin_all <- min(c(X_su_all[,2], X_wi_all[,2])) - pad
ymax_all <- max(c(X_su_all[,2], X_wi_all[,2])) + pad
xlim_all <- c(xmin_all, xmax_all)
ylim_all <- c(ymin_all, ymax_all)

theme_sq <- theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        aspect.ratio = 1)

# 1) Population complète (avant Ω)
df_before <- rbind(
  data.frame(x = X_su_all[,1], y = X_su_all[,2], saison = "été"),
  data.frame(x = X_wi_all[,1], y = X_wi_all[,2], saison = "hiver")
)

p_before <- ggplot(df_before, aes(x, y, color = saison)) +
  geom_point(size = 0.55, alpha = 0.85) +
  scale_color_manual(values = c("été" = "#F28E2B", "hiver" = "#4E79A7")) +
  coord_equal(xlim = xlim_all, ylim = ylim_all, expand = FALSE) +
  theme_sq + labs(title = "AVANT Ω — Été & Hiver", x = "x", y = "y", color = "Saison")

# 2) Fonds Ω (été / hiver)
su_r_df <- as.data.frame(r_su, xy = TRUE); names(su_r_df)[3] <- "omega"
wi_r_df <- as.data.frame(r_wi, xy = TRUE); names(wi_r_df)[3] <- "omega"

p_su_omega <- ggplot(su_r_df, aes(x, y, fill = omega)) +
  geom_raster() +
  stat_contour(aes(z = omega), color = "white", alpha = 0.7,
               breaks = c(0.5, 0.75, 0.9) * max(su_r_df$omega)) +
  scale_fill_viridis_c(name = expression(Omega[été])) +
  coord_equal(xlim = xlim_all, ylim = ylim_all, expand = FALSE) +
  theme_sq + labs(title = "Fond Ω — ÉTÉ", x = "x", y = "y")

p_wi_omega <- ggplot(wi_r_df, aes(x, y, fill = omega)) +
  geom_raster() +
  stat_contour(aes(z = omega), color = "white", alpha = 0.7,
               breaks = c(0.5, 0.75, 0.9) * max(wi_r_df$omega)) +
  scale_fill_viridis_c(name = expression(Omega[hiver])) +
  coord_equal(xlim = xlim_all, ylim = ylim_all, expand = FALSE) +
  theme_sq + labs(title = "Fond Ω — HIVER", x = "x", y = "y")

# 3) Points observés après Ω (avec saison)
df_after <- rbind(
  transform(obs_long, x = x_su, y = y_su, saison = "été")[, c("x","y","saison")],
  transform(obs_long, x = x_wi, y = y_wi, saison = "hiver")[, c("x","y","saison")]
)

p_after <- ggplot(df_after, aes(x, y, color = saison)) +
  geom_point(size = 0.60, alpha = 0.90) +
  scale_color_manual(values = c("été" = "#F28E2B", "hiver" = "#4E79A7")) +
  coord_equal(xlim = xlim_all, ylim = ylim_all, expand = FALSE) +
  theme_sq + labs(title = "APRÈS Ω — Été & Hiver", x = "x", y = "y", color = "Saison")

# 4) Mise en page
(fig_points <- (p_before | p_after) / (p_su_omega | p_wi_omega))


