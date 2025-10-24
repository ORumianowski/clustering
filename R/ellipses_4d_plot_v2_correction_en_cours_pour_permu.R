# ============================================================== # 
# Script 2 — Ellipses combinées (été+hiver) + segments été→hiver 
# Panneaux : VRAI | NAÏF | CORRIGÉ (une seule légende en bas) 
# ============================================================== 

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(mvtnorm)
  library(viridisLite); library(patchwork); library(stringr)
})

# ---------- Helpers ----------

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
    if (any(is.na(S))) {
      Prec_k <- Prec_k + diag(D)*1e-4; S <- solve(Prec_k)
    }
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
  dens <- sapply(1:nrow(mu_mat), function(k) pi_vec[k] * mvtnorm::dmvnorm(X4, mean = mu_mat[k,], sigma = Sigma_list[[k]], log = FALSE))
  max.col(dens, ties.method = "first")
}

ellipse_points <- function(mu, Sigma, r, n = 200){
  ang <- seq(0, 2*pi, length.out = n); circle <- cbind(cos(ang), sin(ang))
  R <- chol(Sigma); pts <- circle %*% t(R) * r; sweep(pts, 2, mu, "+")
}

make_aplats_sets <- function(mu_su, Sig_su, mu_wi, Sig_wi, probs = c(0.50, 0.75, 0.90)){
  r_levels <- sqrt(qchisq(probs, df = 2))
  cov_labels <- sprintf("%d%%", round(100*probs))
  out <- vector("list", 2*nrow(mu_su)*length(r_levels)); idx <- 1L
  for (k in 1:nrow(mu_su)){
    for (i in seq_along(r_levels)){
      pts <- ellipse_points(mu_su[k,], Sig_su[[k]], r_levels[i])
      out[[idx]] <- data.frame(x=pts[,1], y=pts[,2], cluster=factor(k), set="été", level=cov_labels[i]); idx <- idx+1L
    }
    for (i in seq_along(r_levels)){
      pts <- ellipse_points(mu_wi[k,], Sig_wi[[k]], r_levels[i])
      out[[idx]] <- data.frame(x=pts[,1], y=pts[,2], cluster=factor(k), set="hiver", level=cov_labels[i]); idx <- idx+1L
    }
  }
  out <- do.call(rbind, out)
  out$set <- factor(out$set, levels=c("été","hiver"))
  out$level <- factor(out$level, levels=c("50%","75%","90%"))
  out
}

fmt_pi <- function(v) paste0("π = [", paste(sprintf("%.2f", v), collapse=", "), "]")
make_panel_title <- function(lbl) paste0("Ellipses combinées — ", lbl, " (été+hiver)")

# ---------- Fonction de réordonnancement SIMPLIFIÉE ----------

reorder_clusters <- function(mu_mat, Sigma_list, pi_vec) {
  # Trier les clusters par x1 croissant (coordonnée été x1)
  x1_su <- mu_mat[, 1]
  ord <- order(x1_su)
  
  # Réordonner tous les paramètres
  list(
    mu = mu_mat[ord, , drop = FALSE],
    Sigma = Sigma_list[ord],
    pi = pi_vec[ord]
  )
}

# ---------- Fonction principale CORRIGÉE ----------

plot_results_mixture <- function(samples_corr4, samples_naif4, K, X4_all, sel_ete, sel_hiver, X_su_all, X_wi_all, mu_true = NULL, Sigma_true = NULL, pi_true = NULL, z_true = NULL, subplot_widths = c(1,1,1), panel_aspect = 0.9) {
  samples_corr <- as.data.frame(samples_corr4)
  samples_naif <- as.data.frame(samples_naif4)
  
  # Emprise commune
  pad <- 1
  xlim_all <- range(c(X_su_all[,1], X_wi_all[,1])) + c(-pad,pad)
  ylim_all <- range(c(X_su_all[,2], X_wi_all[,2])) + c(-pad,pad)
  
  # Postérieurs
  post_c <- extr_muSigma4D(samples_corr, K, D=4)
  post_n <- extr_muSigma4D(samples_naif, K, D=4)
  pi_c <- get_pi_means(samples_corr, K)
  pi_n <- get_pi_means(samples_naif, K)
  
  # RÉORDONNANCEMENT: Trier par x1 croissant
  ord_c <- reorder_clusters(post_c$mu, post_c$Sigma, pi_c)
  ord_n <- reorder_clusters(post_n$mu, post_n$Sigma, pi_n)
  
  post_c$mu <- ord_c$mu
  post_c$Sigma <- ord_c$Sigma
  pi_c <- ord_c$pi
  
  post_n$mu <- ord_n$mu
  post_n$Sigma <- ord_n$Sigma
  pi_n <- ord_n$pi
  
  # RÉORDONNANCEMENT de la vérité si disponible
  if (!is.null(mu_true) && !is.null(Sigma_true) && !is.null(pi_true)) {
    mu_true_mat <- do.call(rbind, mu_true)
    x1_su_true <- mu_true_mat[, 1]
    ord_true <- order(x1_su_true)
    mu_true <- mu_true[ord_true]
    Sigma_true <- Sigma_true[ord_true]
    pi_true <- pi_true[ord_true]
  }
  
  # Extraire les composantes été/hiver APRÈS réordonnancement
  mu_su_c <- post_c$mu[,1:2,drop=FALSE]; mu_wi_c <- post_c$mu[,3:4,drop=FALSE]
  mu_su_n <- post_n$mu[,1:2,drop=FALSE]; mu_wi_n <- post_n$mu[,3:4,drop=FALSE]
  Sig_su_c <- lapply(post_c$Sigma, \(S) S[1:2,1:2,drop=FALSE])
  Sig_wi_c <- lapply(post_c$Sigma, \(S) S[3:4,3:4,drop=FALSE])
  Sig_su_n <- lapply(post_n$Sigma, \(S) S[1:2,1:2,drop=FALSE])
  Sig_wi_n <- lapply(post_n$Sigma, \(S) S[3:4,3:4,drop=FALSE])
  
  # Segments pour individus observés - UTILISER LES PARAMÈTRES RÉORDONNÉS
  X_obs_any <- X4_all[(sel_ete | sel_hiver), , drop=FALSE]
  
  # Réassigner avec les paramètres réordonnés
  z_corr <- factor(assign_MAP(X_obs_any, post_c$mu, post_c$Sigma, pi_c))
  z_naif <- factor(assign_MAP(X_obs_any, post_n$mu, post_n$Sigma, pi_n))
  
  mk_pairs <- function(X4, zfac) data.frame(lonA=X4[,1], latA=X4[,2], lonB=X4[,3], latB=X4[,4], cluster=zfac)
  
  set.seed(1)
  use_idx <- if (nrow(X_obs_any) > 1000) sample(nrow(X_obs_any), 1000) else seq_len(nrow(X_obs_any))
  
  # Créer les segments avec les assignments réordonnés
  pairs_corr <- mk_pairs(X_obs_any[use_idx,], z_corr[use_idx])
  pairs_naif <- mk_pairs(X_obs_any[use_idx,], z_naif[use_idx])
  
  # Ellipses avec les paramètres réordonnés
  polys_corr <- make_aplats_sets(mu_su_c, Sig_su_c, mu_wi_c, Sig_wi_c)
  polys_naif <- make_aplats_sets(mu_su_n, Sig_su_n, mu_wi_n, Sig_wi_n)
  
  # Palettes (cohérentes grâce au réordonnancement)
  cluster_palette <- setNames(rocket(K, begin=0.10, end=0.75, direction=-1), as.character(1:K))
  coverage_alphas <- c("50%"=.50, "75%"=.20, "90%"=.10)
  
  # Fabrique d'un panneau
  build_pairs_plot <- function(pairs_df, polys_df, title, subtitle){
    ggplot() +
      geom_segment(data = pairs_df, aes(x=lonA, y=latA, xend=lonB, yend=latB, color=cluster), 
                   alpha=.06, linewidth=.5, show.legend = FALSE) +
      geom_polygon(data = polys_df, aes(x=x, y=y, group=interaction(cluster, level, set), 
                                        fill=cluster, alpha=level), color=NA) +
      scale_fill_manual(values=cluster_palette, name="Clusters", 
                        guide = guide_legend(override.aes = list(alpha = 1))) +
      scale_color_manual(values=cluster_palette, guide = "none") +
      scale_alpha_manual(values=coverage_alphas, name="Contours prob.") +
      coord_fixed(ratio = panel_aspect, xlim=xlim_all, ylim=ylim_all, expand=FALSE) +
      theme_minimal(base_size=12) +
      theme(plot.title = element_text(face="bold")) +
      labs(title=title, subtitle=subtitle, x="x", y="y")
  }
  
  # Panneaux NAÏF & CORRIGÉ
  p_naif <- build_pairs_plot(pairs_naif, polys_naif, make_panel_title("NAÏF"), fmt_pi(pi_n))
  p_corr <- build_pairs_plot(pairs_corr, polys_corr, make_panel_title("CORRIGÉ"), fmt_pi(pi_c))
  
  # Panneau VRAI (si mu/Sigma fournis)
  has_truth <- !is.null(mu_true) && !is.null(Sigma_true)
  if (has_truth) {
    mu_true_mat <- do.call(rbind, mu_true)
    mu_su_true <- mu_true_mat[,1:2,drop=FALSE]
    mu_wi_true <- mu_true_mat[,3:4,drop=FALSE]
    Sig_su_true <- lapply(Sigma_true, \(S) S[1:2,1:2,drop=FALSE])
    Sig_wi_true <- lapply(Sigma_true, \(S) S[3:4,3:4,drop=FALSE])
    polys_true <- make_aplats_sets(mu_su_true, Sig_su_true, mu_wi_true, Sig_wi_true)
    
    # Préparer les segments pour VRAI
    if (!is.null(z_true)) {
      # Réordonner z_true si nécessaire
      if (!is.null(pi_true)) {
        # On suppose que z_true est déjà dans le bon ordre après réordonnancement des paramètres
        pairs_true_df <- data.frame(
          lonA = X_obs_any[use_idx, 1], latA = X_obs_any[use_idx, 2],
          lonB = X_obs_any[use_idx, 3], latB = X_obs_any[use_idx, 4],
          cluster = factor(z_true[use_idx])
        )
      } else {
        # Si pas de pi_true, utiliser les assignments MAP avec les paramètres vrais réordonnés
        z_true_map <- factor(assign_MAP(X_obs_any, mu_true_mat, Sigma_true, pi_true))
        pairs_true_df <- data.frame(
          lonA = X_obs_any[use_idx, 1], latA = X_obs_any[use_idx, 2],
          lonB = X_obs_any[use_idx, 3], latB = X_obs_any[use_idx, 4],
          cluster = z_true_map[use_idx]
        )
      }
    } else {
      # Si z_true pas disponible, utiliser les assignments MAP
      z_true_map <- factor(assign_MAP(X_obs_any, mu_true_mat, Sigma_true, pi_true))
      pairs_true_df <- data.frame(
        lonA = X_obs_any[use_idx, 1], latA = X_obs_any[use_idx, 2],
        lonB = X_obs_any[use_idx, 3], latB = X_obs_any[use_idx, 4],
        cluster = z_true_map[use_idx]
      )
    }
    
    # Sous-titre π (vraies valeurs si dispo)
    pi_true_sub <- NULL
    if (!is.null(pi_true)) {
      pi_true_sub <- as.numeric(pi_true)
    } else if (!is.null(z_true)) {
      tabz <- prop.table(table(z_true))
      pi_true_sub <- replace(numeric(K), 1:K, as.numeric(tabz[as.character(1:K)]))
      pi_true_sub[is.na(pi_true_sub)] <- 0
    }
    
    p_true <- build_pairs_plot(
      pairs_true_df, polys_true, make_panel_title("VRAI"), 
      if (!is.null(pi_true_sub)) fmt_pi(pi_true_sub) else "π (simulation)"
    )
    
    out <- (p_true | p_naif | p_corr) + plot_layout(guides = "collect", widths = subplot_widths) & 
      theme(legend.position = "bottom")
  } else {
    out <- (p_naif | p_corr) + plot_layout(guides = "collect", widths = subplot_widths[1:2]) & 
      theme(legend.position = "bottom")
  }
  out
}

# ------------------ EXÉCUTION ------------------

# Ce bloc suppose que les objets suivants existent déjà dans ton environnement :
# samples_corr4, samples_naif4, K,
# X4_all, sel_ete, sel_hiver, X_su_all, X_wi_all,
# (optionnels) mu_true, Sigma_true, pi_true, z

subplot_widths <- c(1,1,1) # largeurs relatives des trois panneaux
panel_aspect <- 0.9 # ratio hauteur/largeur de chaque panneau

res_plots <- plot_results_mixture(
  samples_corr4, samples_naif4, K,
  X4_all, sel_ete, sel_hiver, X_su_all, X_wi_all,
  # passe ces arguments s'ils existent chez toi :
  mu_true = if (exists("mu_true")) mu_true else NULL,
  Sigma_true = if (exists("Sigma_true")) Sigma_true else NULL,
  pi_true = if (exists("pi_true")) pi_true else NULL,
  z_true = if (exists("z")) z else NULL,
  subplot_widths = subplot_widths,
  panel_aspect = panel_aspect
)

print(res_plots)