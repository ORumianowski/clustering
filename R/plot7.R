############################################################
# Figure publication-ready : GMM ellipses + CMR/GPS
# - CMR : cercles blancs bordés de noir
# - GPS : losanges blancs bordés de noir
# - Ellipses : couleurs par cluster (viridis), très translucides
# - Centroids : croix colorées, halo blanc
# - Légendes séparées (Data source / Clusters / Probability contours)
# - Ajout : filtrage d’un cluster donné (k_sel)
############################################################

# ============== 0) Packages ==============
library(nimble)
library(coda)
library(ggplot2)
library(dplyr)
library(viridisLite)
library(patchwork)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)

set.seed(123)

# ============== 1) Données ==============
nb_clusters <- 3

load("dataset/soft_connectivity.data.RData")
gps <- read.csv("dataset/curlew_metadata_22-02-2024_coordinates.csv",
                sep = ",", dec = ".", header = TRUE)

X_cmr <- cbind(
  soft_connectivity.data$longitude_breed,
  soft_connectivity.data$latitude_breed,
  soft_connectivity.data$longitude_wint,
  soft_connectivity.data$latitude_wint
) |> as.data.frame() |> na.omit() |> mutate(source = "cmr")

X_gps <- cbind(
  gps$Breeding.site.longitude,
  gps$Breeding.site.latitude,
  gps$Wintering.site.1.longitude,
  gps$Wintering.site.1.latitude
) |> as.data.frame() |> na.omit() |> mutate(source = "gps")

X_full <- rbind(X_cmr, X_gps)
X      <- X_full[, 1:4]
colnames(X) <- c("lon_breed","lat_breed","lon_wint","lat_wint")

# ============== 2) Modèle NIMBLE ==============
code_connectivity_only <- nimbleCode({
  pi[1:K] ~ ddirch(alpha[1:K])
  for(k in 1:K){
    for(j in 1:d){
      mu[k,j]    ~ dnorm(0, var = 10)
      sigma[k,j] ~ dinvgamma(shape = 2.1, scale = 1.1)
      prec[k,j] <- 1 / sigma[k,j]
    }
  }
  for(i in 1:N){
    cl[i] ~ dcat(pi[1:K])
    for(j in 1:d){
      X[i,j] ~ dnorm(mu[cl[i],j], prec[cl[i],j])
    }
  }
})

N <- nrow(X); d <- ncol(X); K <- nb_clusters
d_means  <- colMeans(X, na.rm = TRUE)
init_mu  <- t(replicate(K, d_means + rnorm(d, 0, 1)))

constants      <- list(K = K, N = N, d = d, alpha = rep(1, K))
data_conn_only <- list(X = as.matrix(X))
inits_conn     <- list(
  cl = sample(1:K, N, replace = TRUE),
  mu = init_mu,
  sigma = matrix(1, nrow = K, ncol = d),
  pi = rep(1/K, K)
)

run_nimble <- function(code, constants, data, inits, monitors,
                       niter = 6000, nburnin = 2000, thin = 4,
                       progressBar = TRUE) {
  model  <- nimbleModel(code, constants = constants, data = data, inits = inits, calculate = TRUE)
  cmodel <- compileNimble(model)
  conf   <- configureMCMC(model, monitors = monitors)
  mcmc   <- buildMCMC(conf)
  cmcmc  <- compileNimble(mcmc, project = model)
  samp   <- runMCMC(cmcmc, niter = niter, nburnin = nburnin, thin = thin, progressBar = progressBar)
  list(model = model, cmodel = cmodel, samples = as.mcmc(samp))
}

cat("Running connectivity-only model...\n")
res_conn <- run_nimble(code_connectivity_only, constants, data_conn_only, inits_conn,
                       monitors = c("pi","mu","sigma","cl"),
                       niter = 6000, nburnin = 2000, thin = 4)

# ============== 3) Récupération des paramètres (posterior means) ==============
samples_conn_df <- as.data.frame(as.matrix(res_conn$samples))

# Moyennes des mu
mu_cols <- grep("^mu\\[", colnames(samples_conn_df), value = TRUE)
k_idx   <- as.integer(sub("^mu\\[(\\d+),\\s*(\\d+)\\]$", "\\1", mu_cols))
d_idx   <- as.integer(sub("^mu\\[(\\d+),\\s*(\\d+)\\]$", "\\2", mu_cols))
K <- max(k_idx); D <- max(d_idx)

mu_post <- matrix(NA_real_, nrow = K, ncol = D)
for (c in seq_along(mu_cols)) {
  mu_post[k_idx[c], d_idx[c]] <- mean(samples_conn_df[[mu_cols[c]]])
}
colnames(mu_post) <- c("lon_breed","lat_breed","lon_wint","lat_wint")

Sigma_post <- vector("list", K)
for (k in 1:K) {
  sigmas <- numeric(D)
  for (jj in 1:D) {
    colname <- paste0("sigma[", k, ", ", jj, "]")
    sigmas[jj] <- mean(samples_conn_df[[colname]])
  }
  Sigma_post[[k]] <- diag(sigmas, D, D)
}

# Attribution MAP des individus
cl_cols <- grep("^cl\\[", colnames(samples_conn_df), value = TRUE)
cl_mat <- samples_conn_df[, cl_cols]
cl_map <- apply(cl_mat, 2, function(v) as.integer(names(sort(table(v), decreasing = TRUE)[1])))
X_full$cluster <- cl_map

# ============== 4) Préparation des données pour affichage ==============
# Palette clusters
cluster_levels  <- factor(seq_len(K), levels = as.character(seq_len(K)))
cluster_palette <- setNames(
  viridisLite::rocket(K, begin = 0.10, end = 0.75, direction = -1),
  levels(cluster_levels)
)

# Ellipses
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
        cluster  = factor(k, levels = levels(cluster_levels)),
        coverage = factor(cov_labels[i], 
                          levels = c("50%","75%","90%","95%","99%"))
      )
      idx <- idx + 1L
    }
  }
  dplyr::bind_rows(out)
}

mu_summer    <- mu_post[,1:2, drop = FALSE]
mu_winter    <- mu_post[,3:4, drop = FALSE]
Sigma_summer <- lapply(Sigma_post, \(M) M[1:2,1:2])
Sigma_winter <- lapply(Sigma_post, \(M) M[3:4,3:4])

polys_summer <- make_aplats(split(mu_summer, row(mu_summer)), Sigma_summer)
polys_winter <- make_aplats(split(mu_winter, row(mu_winter)), Sigma_winter)

# ============== 5) Thème & base map ==============
world <- ne_countries(scale = "medium", returnclass = "sf")
base_map <- list(
  geom_sf(data = world, fill = "grey95", color = "grey80", linewidth = 0.3)
)

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

coverage_alphas <- c("50%" = 0.20,
                     "75%" = 0.12,
                     "90%" = 0.08,
                     "95%" = 0.05,
                     "99%" = 0.03)

# ============== 6) Préparer un cluster choisi sur une seule carte ==============
k_sel <- 2  # <-- changer ici pour le cluster désiré

# Couples breeding/wintering
df_pairs <- X_full |> 
  filter(cluster == k_sel) |> 
  mutate(id = row_number()) |> 
  select(id, source, lon_breed, lat_breed, lon_wint, lat_wint)

# Points
df_points <- df_pairs |>
  tidyr::pivot_longer(
    cols = c(lon_breed, lat_breed, lon_wint, lat_wint),
    names_to = c(".value","site"),
    names_pattern = "(lon|lat)_(.*)"
  ) |>
  mutate(site = ifelse(site == "breed", "Breeding", "Wintering"))

# Segments
df_segments <- df_pairs |> 
  transmute(id, lon_b_
            
