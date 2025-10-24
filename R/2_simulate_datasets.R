# =============================================================
# Script 2/4 — Simulation des jeux de données (3 contextes) + plots
# =============================================================
# IMPORTANT: La simulation utilise TOUJOURS Ω "original" (été & hiver).
#            Les versions dégradées d'Ω (deg1/deg2) ne sont utilisées
#            QUE dans l'analyse (Script 3) pour tester l'impact d'une
#            carte d'effort imparfaite. Ici on produit 3 datasets (C1,C2,C3).
# =============================================================

suppressPackageStartupMessages({
  library(MASS)
  library(terra)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

# --- IMPORTANT: réutiliser Script 1 ---
if (!exists("get_context_params")) {
  source("R/1_contexts_omegas.R")
}

# ---------- Simulation 4D soumise à Ω ORIGINAL uniquement ----------
simulate_dataset <- function(ctx=c("C1","C2","C3"),
                             N_tot=5000, pi_true=c(0.20,0.35,0.45), seed=123, n_pix_side=60){
  stopifnot(length(pi_true) == 3)
  set.seed(seed)
  ctx <- match.arg(ctx)
  pars <- get_context_params(ctx)
  mu_true <- pars$mu_true; Sigma_true <- pars$Sigma_true; K <- length(mu_true); D <- 4
  
  # Attribution des clusters
  z <- sample(1:K, N_tot, replace = TRUE, prob = pi_true)
  
  # Échantillonnage 4D
  X4_all <- t(vapply(1:N_tot, function(i) MASS::mvrnorm(1, mu = mu_true[[z[i]]], Sigma = Sigma_true[[z[i]]]), numeric(D)))
  colnames(X4_all) <- c("x_su","y_su","x_wi","y_wi")
  X_su_all <- X4_all[,1:2, drop = FALSE]
  X_wi_all <- X4_all[,3:4, drop = FALSE]
  
  # Ω ORIGINAL — on fixe la version utilisée pour la simulation
  om <- make_omega_pair("original", n_pix_side = n_pix_side)
  r_su <- om$summer; r_wi <- om$winter
  
  Omega_su_all <- terra::extract(r_su, X_su_all)[,1]
  Omega_wi_all <- terra::extract(r_wi, X_wi_all)[,1]
  Omega_su_all[is.na(Omega_su_all)] <- 0
  Omega_wi_all[is.na(Omega_wi_all)] <- 0
  
  sel_ete   <- runif(N_tot) < Omega_su_all
  sel_hiver <- runif(N_tot) < Omega_wi_all
  
  obs_ete <- data.frame(
    id = which(sel_ete), saison = "ete",
    x_su = X_su_all[sel_ete,1], y_su = X_su_all[sel_ete,2],
    x_wi = X_wi_all[sel_ete,1], y_wi = X_wi_all[sel_ete,2],
    Omega = Omega_su_all[sel_ete], z = z[sel_ete]
  )
  obs_hiver <- data.frame(
    id = which(sel_hiver), saison = "hiver",
    x_su = X_su_all[sel_hiver,1], y_su = X_su_all[sel_hiver,2],
    x_wi = X_wi_all[sel_hiver,1], y_wi = X_wi_all[sel_hiver,2],
    Omega = Omega_wi_all[sel_hiver], z = z[sel_hiver]
  )
  
  obs_long <- rbind(obs_ete, obs_hiver)
  obs_long$ctx <- ctx
  obs_long$omega_used_for_sim <- "original"
  
  list(
    data_long = obs_long,
    X4_all = X4_all,
    z_all = z,
    mu_true = mu_true,
    Sigma_true = Sigma_true,
    r_su = r_su,
    r_wi = r_wi
  )
}

# ---------- Plot des nuages de points (par saison), couleur = cluster vrai ----------
plot_dataset_points <- function(sim, max_points = 5000){
  xmin <- if (exists(".XMIN_ALL")) .XMIN_ALL else min(c(sim$data_long$x_su, sim$data_long$x_wi), na.rm=TRUE) - 5
  xmax <- if (exists(".XMAX_ALL")) .XMAX_ALL else max(c(sim$data_long$x_su, sim$data_long$x_wi), na.rm=TRUE) + 5
  ymin <- if (exists(".YMIN_ALL")) .YMIN_ALL else min(c(sim$data_long$y_su, sim$data_long$y_wi), na.rm=TRUE) - 5
  ymax <- if (exists(".YMAX_ALL")) .YMAX_ALL else max(c(sim$data_long$y_su, sim$data_long$y_wi), na.rm=TRUE) + 5
  
  df <- sim$data_long
  if (nrow(df) > max_points) df <- df[sample.int(nrow(df), max_points), ]
  df$cluster <- factor(df$z)
  
  df_su <- df |> dplyr::mutate(x = x_su, y = y_su) |> dplyr::filter(saison == "ete")
  df_wi <- df |> dplyr::mutate(x = x_wi, y = y_wi) |> dplyr::filter(saison == "hiver")
  
  p_su <- ggplot(df_su, aes(x, y, color=cluster)) + geom_point(alpha=.6, size=0.8) +
    coord_fixed(xlim=c(xmin,xmax), ylim=c(ymin,ymax), expand=FALSE) +
    theme_minimal() + labs(title=sprintf("%s — Été (Ω sim: original)", unique(df$ctx)), x="x", y="y")
  p_wi <- ggplot(df_wi, aes(x, y, color=cluster)) + geom_point(alpha=.6, size=0.8) +
    coord_fixed(xlim=c(xmin,xmax), ylim=c(ymin,ymax), expand=FALSE) +
    theme_minimal() + labs(title=sprintf("%s — Hiver (Ω sim: original)", unique(df$ctx)), x="x", y="y")
  
  p_su / p_wi
}

# ---------- Génération des 3 datasets (C1,C2,C3) ----------
make_all_datasets_contexts <- function(N_tot=5000, seed_base=123, n_pix_side=60){
  ctxs <- c("C1","C2","C3")
  sims <- vector("list", length(ctxs))
  for (i in seq_along(ctxs)){
    sims[[i]] <- simulate_dataset(ctxs[i], N_tot=N_tot, seed = seed_base + i, n_pix_side = n_pix_side)
  }
  names(sims) <- ctxs
  sims
}

# ---------- Visualisation multi-contextes ----------
plot_all_datasets <- function(sims){
  plots <- lapply(sims, plot_dataset_points)
  wrap_plots(plots, ncol = 3)
}

# =============================================================
# Exemple d'utilisation :
#   Génère les 3 jeux (C1,C2,C3) avec Ω original et affiche les figures.
# =============================================================

sims_ctx <- make_all_datasets_contexts(N_tot = 5000, seed_base = 123, n_pix_side = 60)
p <- plot_all_datasets(sims_ctx)
print(p)

# ---- FIN Script 2 ----