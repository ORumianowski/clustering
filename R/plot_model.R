
library(ggplot2)
library(ellipse)
library(gridExtra)

# Fonction pour créer un plot 2D pour deux dimensions données
plot_2D_clusters <- function(d1, d2, X, z_true, mu_est_list, Sigma_est_list, true_means, true_covs, K) {
  # Points
  df_points <- data.frame(
    x = X[, d1],
    y = X[, d2],
    cluster = factor(z_true)
  )
  
  # Ellipses estimées
  ellipses_est <- do.call(rbind, lapply(1:K, function(k) {
    e <- ellipse(Sigma_est_list[[k]][c(d1,d2), c(d1,d2)],
                 centre = mu_est_list[[k]][c(d1,d2)],
                 level = 0.95)
    data.frame(e, cluster = factor(k), type = "Estimated")
  }))
  
  # Ellipses vraies
  ellipses_true <- do.call(rbind, lapply(1:K, function(k) {
    e <- ellipse(true_covs[[k]][c(d1,d2), c(d1,d2)],
                 centre = true_means[[k]][c(d1,d2)],
                 level = 0.95)
    data.frame(e, cluster = factor(k), type = "True")
  }))
  
  ellipses_all <- rbind(ellipses_est, ellipses_true)
  
  # Moyennes estimées
  mu_points <- do.call(rbind, lapply(1:K, function(k) {
    data.frame(x = mu_est_list[[k]][d1], y = mu_est_list[[k]][d2],
               cluster = factor(k))
  }))
  
  # Plot
  p <- ggplot() +
    geom_point(data = df_points, aes(x = x, y = y, color = cluster), alpha = 0.6) +
    geom_path(data = ellipses_all, aes(x = x, y = y, color = cluster, linetype = type), size = 1) +
    geom_point(data = mu_points, aes(x = x, y = y, color = cluster), shape = 4, size = 4, stroke = 2) +
    scale_linetype_manual(values = c("Estimated" = "solid", "True" = "dashed")) +
    labs(x = paste("Dimension", d1), y = paste("Dimension", d2)) +
    theme_minimal() +
    theme(legend.position = "right")
  
  return(p)
}


# --- Créer les listes des moyennes et covariances estimées ---
mu_est_list <- list()
Sigma_est_list <- list()
for (k in 1:K) {
  # Moyenne estimée du cluster k
  mu_idx <- grep(paste0("mu\\[", k, ","), colnames(samples))
  mu_est <- colMeans(samples[, mu_idx])
  mu_est_list[[k]] <- mu_est
  
  # Précision estimée et conversion en covariance
  Prec_idx <- grep(paste0("Prec\\[", k, ","), colnames(samples))
  Prec_est <- matrix(colMeans(samples[, Prec_idx]), d, d)
  Sigma_est_list[[k]] <- solve(Prec_est)  # covariance estimée
}


# ============================================================
# Créer les deux plots
# ============================================================
p1 <- plot_2D_clusters(1, 2, X, z_true, mu_est_list, Sigma_est_list, true_means, true_covs, K)
p2 <- plot_2D_clusters(3, 4, X, z_true, mu_est_list, Sigma_est_list, true_means, true_covs, K)

# Affichage côte à côte
grid.arrange(p1, p2, ncol = 2, top = "Projection 2D des clusters (Dimensions 1-2 et 3-4)")
