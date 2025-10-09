library(MASS)
library(dplyr)
library(ggplot2)
library(viridisLite)

set.seed(123)

# --- 1) Simulation 2D
K <- 3
mu_list <- list(c(5,50), c(10,52), c(7,54))
Sigma_list <- list(
  matrix(c(3, 1, 1, 2), 2, 2),
  matrix(c(2, -0.5, -0.5, 2), 2, 2),
  matrix(c(2, 0.8, 0.8, 3), 2, 2)
)
n_per_cluster <- 50

points <- do.call(rbind, lapply(1:K, function(k) {
  pts <- MASS::mvrnorm(n_per_cluster, mu = mu_list[[k]], Sigma = Sigma_list[[k]])
  data.frame(
    lon = pts[,1], lat = pts[,2],
    cluster = factor(k, levels = as.character(1:K)),
    source  = sample(c("cmr","gps"), n_per_cluster, replace = TRUE)
  )
}))

# --- 2) Projection PCA -> 1D
X <- as.matrix(points[,c("lon","lat")])
pca <- prcomp(X, center = TRUE, scale. = FALSE)
v1  <- pca$rotation[,1]
points$proj <- as.numeric(X %*% v1)

# --- 3) Paramètres 1D par cluster (fit)
stats_df <- points %>%
  group_by(cluster) %>%
  summarise(mu = mean(proj), sd = sd(proj), .groups = "drop")

# --- 4) Courbes normales théoriques
x_grid <- seq(min(points$proj) - 1, max(points$proj) + 1, length.out = 500)
dens_df <- do.call(rbind, lapply(1:K, function(k) {
  mu_k <- stats_df$mu[k]
  sd_k <- stats_df$sd[k]
  data.frame(
    x = x_grid,
    y = dnorm(x_grid, mean = mu_k, sd = sd_k),
    cluster = factor(k, levels = levels(points$cluster))
  )
}))

# --- 5) Bandes de quantiles sous les normales
coverage_probs  <- c(0.50, 0.70, 0.90, 0.95)
coverage_labels <- paste0(100*coverage_probs, "%")
coverage_alphas <- c("50%" = 0.20,
                     "70%" = 0.14,
                     "90%" = 0.10,
                     "95%" = 0.06)

build_band_polygon <- function(mu, sd, cov_prob, cl_label) {
  z <- qnorm(1 - (1 - cov_prob)/2)
  q_low  <- mu - z*sd
  q_high <- mu + z*sd
  x_all  <- seq(q_low, q_high, length.out = 200)
  y_all  <- dnorm(x_all, mean = mu, sd = sd)
  data.frame(
    x = c(x_all, rev(x_all)),
    y = c(y_all, rep(0, length(x_all))),
    cluster  = cl_label,
    coverage = factor(paste0(100*cov_prob, "%"), levels = coverage_labels)
  )
}

bands_list <- list()
idx <- 1L
for (k in 1:K) {
  mu_k <- stats_df$mu[k]
  sd_k <- stats_df$sd[k]
  cl   <- stats_df$cluster[k]
  for (cov_prob in coverage_probs) {
    bands_list[[idx]] <- build_band_polygon(mu_k, sd_k, cov_prob, cl)
    idx <- idx + 1L
  }
}
bands_df <- do.call(rbind, bands_list)

# --- 6) Palette
cluster_palette <- setNames(
  viridisLite::rocket(K, begin = 0.10, end = 0.75, direction = -1),
  levels(points$cluster)
)

# --- 7) Plot final
max_y <- max(dens_df$y)
p1d <- ggplot() +
  # Bandes de quantiles
  geom_polygon(data = bands_df,
               aes(x = x, y = y, fill = cluster, alpha = coverage,
                   group = interaction(cluster, coverage)),
               color = NA) +
  # Courbes normales
  geom_line(data = dens_df,
            aes(x = x, y = y, color = cluster), linewidth = 1) +
  # Points projetés
  geom_jitter(data = points,
              aes(x = proj, y = -0.06*max_y, color = cluster, fill = cluster, shape = source),
              height = 0.02*max_y, size = 1.9, alpha = 0.8, stroke = 0.7) +
  coord_cartesian(ylim = c(-0.12*max_y, 1.05*max_y), expand = FALSE) +
  scale_fill_manual(values = cluster_palette, name = "Cluster") +
  scale_color_manual(values = cluster_palette, name = "Cluster") +
  scale_alpha_manual(values = coverage_alphas, name = "Quantile bands") +
  scale_shape_manual(values = c("cmr" = 21, "gps" = 23), name = "Data source") +
  labs(title = "Projection 1D : Normales ajustées + bandes de quantiles",
       x = "Coordonnée projetée (axe PCA 1)", y = "Densité") +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.box      = "vertical",
    axis.title      = element_text(face = "bold")
  )

print(p1d)
