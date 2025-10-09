library(MASS)
library(dplyr)
library(ggplot2)
library(viridisLite)

set.seed(123)

# ============== 1) Simulation: couples (A/B) par cluster ==============
K <- 3
n_per_cluster <- 30

# Moyennes des groupes en haut (A)
mu_list <- list(
  c(3, 9),   # cluster 1
  c(6, 9),   # cluster 2
  c(2.1, 9)    # cluster 3
)

# Covariances serrées (peu de dispersion)
Sigma_list <- list(
  matrix(c(0.25,  0.05, 0.05, 0.25), 2, 2),
  matrix(c(0.25, -0.05,-0.05, 0.25), 2, 2),
  matrix(c(0.25,  0.05, 0.05, 0.25), 2, 2)
)

# Moyennes en bas (B)
# -> clusters 1 et 2 : proches mais pas identiques (recouvrement fort)
# -> cluster 3 : plus à l’écart
mu_list_B <- list(
  c(5, -2),   # cluster 1
  c(5.8, -2.2), # cluster 2 (proche de 1, overlap)
  c(2, -4)    # cluster 3 (séparé)
)

# Génération des couples A/B
pairs <- do.call(rbind, lapply(1:K, function(k) {
  A <- MASS::mvrnorm(n_per_cluster, mu = mu_list[[k]],   Sigma = Sigma_list[[k]])
  B <- MASS::mvrnorm(n_per_cluster, mu = mu_list_B[[k]], Sigma = Sigma_list[[k]])
  data.frame(
    id      = paste0("cl", k, "_", seq_len(n_per_cluster)),
    cluster = factor(k, levels = as.character(1:K)),
    lonA = A[,1], latA = A[,2],
    lonB = B[,1], latB = B[,2]
  )
}))

# ============== 2) Ellipses gaussiennes (A et B) ==============
ellipse_points <- function(mu, Sigma, r, n = 200) {
  ang <- seq(0, 2*pi, length.out = n)
  circle <- cbind(cos(ang), sin(ang))
  R <- chol(Sigma)
  pts <- circle %*% t(R) * r
  sweep(pts, 2, mu, "+")
}
make_aplats <- function(mu_list, Sigma_list, probs = c(0.50, 0.75, 0.90), set_label = "A") {
  r_levels   <- sqrt(qchisq(probs, df = 2))
  cov_labels <- sprintf("%d%%", round(100*probs))
  out <- vector("list", length(mu_list) * length(r_levels))
  idx <- 1L
  for (k in seq_along(mu_list)) {
    for (i in seq_along(r_levels)) {
      pts <- ellipse_points(mu_list[[k]], Sigma_list[[k]], r_levels[i])
      out[[idx]] <- data.frame(
        x = pts[,1], y = pts[,2],
        cluster  = factor(k, levels = as.character(1:K)),
        coverage = factor(cov_labels[i], levels = c("50%","75%","90%")),
        set = factor(set_label, levels = c("A","B"))
      )
      idx <- idx + 1L
    }
  }
  dplyr::bind_rows(out)
}

polys_A <- make_aplats(mu_list,   Sigma_list, set_label = "A")
polys_B <- make_aplats(mu_list_B, Sigma_list, set_label = "B")
polys   <- dplyr::bind_rows(polys_A, polys_B)

# ============== 3) Palette & transparences ==============
cluster_palette <- setNames(
  viridisLite::rocket(K, begin = 0.10, end = 0.75, direction = -1),
  levels(pairs$cluster)
)
coverage_alphas <- c("50%" = 0.20, "75%" = 0.14, "90%" = 0.08)

# ============== 4) Thème ==============
theme_pub <- theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.title       = element_text(face = "bold"),
    axis.text        = element_text(color = "grey10"),
    plot.title       = element_text(face = "bold", hjust = 0.5),
    legend.position  = "bottom",
    legend.box       = "vertical",
    legend.title     = element_text(face = "bold")
  )

# ============== 5) Tracé (segments + ellipses) ==============
build_map <- function(pairs_df, polys_df, title) {
  ggplot() +
    # Segments A -> B
    geom_segment(data = pairs_df,
                 aes(x = lonA, y = latA, xend = lonB, yend = latB, color = cluster),
                 alpha = 0.35, linewidth = 0.7) +
    # Ellipses A & B
    geom_polygon(data = polys_df,
                 aes(x = x, y = y,
                     group = interaction(cluster, coverage, set),
                     fill  = cluster, alpha = coverage),
                 color = NA) +
    scale_fill_manual(values = cluster_palette, name = "Clusters") +
    scale_color_manual(values = cluster_palette, name = "Clusters") +
    scale_alpha_manual(values = coverage_alphas, name = "Probability contours") +
    guides(
      color = guide_legend(order = 1),
      fill  = guide_legend(order = 2),
      alpha = guide_legend(order = 3)
    ) +
    labs(title = title, x = "Longitude", y = "Latitude") +
    theme_pub
}

# ============== 6) Affichage ==============
p <- build_map(
  pairs, polys,
  "Bas : clusters 1 & 2 fortement recouvrants ; haut : distant"
)
print(p)
