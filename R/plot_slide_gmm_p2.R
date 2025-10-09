library(MASS)        
library(dplyr)
library(ggplot2)
library(viridisLite)

set.seed(123)

# ============== 1) Simulation des données ==============
K <- 3  
D <- 2  

mu_list <- list(
  c(5, 50),
  c(10, 52),
  c(7, 54)
)

Sigma_list <- list(
  matrix(c(3, 1, 1, 2), 2, 2),
  matrix(c(2, -0.5, -0.5, 2), 2, 2),
  matrix(c(2, 0.8, 0.8, 3), 2, 2)
)

n_per_cluster <- 50
points <- do.call(rbind, lapply(1:K, function(k) {
  pts <- MASS::mvrnorm(n_per_cluster, mu = mu_list[[k]], Sigma = Sigma_list[[k]])
  data.frame(lon = pts[,1], lat = pts[,2],
             cluster = factor(k),
             source = sample(c("cmr","gps"), n_per_cluster, replace = TRUE))
}))

# ============== 2) Ellipses gaussiennes ==============
ellipse_points <- function(mu, Sigma, r, n = 200) {
  ang <- seq(0, 2*pi, length.out = n)
  circle <- cbind(cos(ang), sin(ang))
  R <- chol(Sigma)
  pts <- circle %*% t(R) * r
  sweep(pts, 2, mu, FUN = "+")
}

make_aplats <- function(mu_list, Sigma_list,
                        probs = c(0.50, 0.75, 0.90, 0.95, 0.99)) {
  r_levels   <- sqrt(qchisq(probs, df = 2))
  cov_labels <- sprintf("%d%%", round(100*probs))
  out <- vector("list", length(mu_list) * length(r_levels))
  idx <- 1L
  for (k in seq_along(mu_list)) {
    for (i in seq_along(r_levels)) {
      pts <- ellipse_points(mu_list[[k]], Sigma_list[[k]], r_levels[i])
      out[[idx]] <- data.frame(
        x        = pts[,1],
        y        = pts[,2],
        cluster  = factor(k, levels = as.character(1:K)),
        coverage = factor(cov_labels[i],
                          levels = c("50%","75%","90%","95%","99%"))
      )
      idx <- idx + 1L
    }
  }
  dplyr::bind_rows(out)
}

polys <- make_aplats(mu_list, Sigma_list)

# Centroides
df_cents <- data.frame(
  x = sapply(mu_list, `[`, 1),
  y = sapply(mu_list, `[`, 2),
  cluster = factor(1:K)
)

# Palette
cluster_palette <- setNames(
  viridisLite::rocket(K, begin = 0.10, end = 0.75, direction = -1),
  as.character(1:K)
)

coverage_alphas <- c("50%" = 0.20,
                     "75%" = 0.12,
                     "90%" = 0.08,
                     "95%" = 0.05,
                     "99%" = 0.03)

# ============== 3) Thème simplifié ==============
theme_pub <- theme_minimal(base_size = 12) +
  theme(
    text             = element_text(family = "Helvetica"),
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

# ============== 4) Fonction de tracé ==============
build_map <- function(points_df, polys_df, cents_df, title) {
  ggplot() +
    
    # Points colorés par cluster
    geom_point(data = points_df,
               aes(x = lon, y = lat, shape = source, color = cluster),
               size = 2.3, stroke = 0.8, alpha = 0.7) +
    
    # Ellipses
    geom_polygon(data = polys_df,
                 aes(x = x, y = y,
                     group = interaction(cluster, coverage),
                     fill  = cluster,
                     alpha = coverage), 
                 color = NA) +
    scale_fill_manual(values = cluster_palette, name = "Clusters") +
    scale_color_manual(values = cluster_palette, name = "Clusters") +
    scale_alpha_manual(values = coverage_alphas, name = "Probability contours") +
    
    # Centroides
    geom_point(data = cents_df,
               aes(x = x, y = y, color = cluster),
               shape = 3, size = 2.8, stroke = 1.2, show.legend = FALSE) +
    
    # Légende sources
    scale_shape_manual(values = c("cmr" = 21, "gps" = 23), name = "Data source") +
    guides(
      shape = guide_legend(order = 1,
                           override.aes = list(size = 3, stroke = 0.8)),
      color = guide_legend(order = 2),
      fill  = guide_legend(order = 2),
      alpha = guide_legend(order = 3)
    ) +
    
    labs(title = title, x = "Longitude", y = "Latitude") +
    theme_pub
}

# ============== 5) Affichage ==============
p <- build_map(points, polys, df_cents, "Simulation: 3 groupes recoupés (points colorés par cluster)")
print(p)
