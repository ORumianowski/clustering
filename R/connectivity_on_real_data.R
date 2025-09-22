
# Script


library(nimble)
library(MASS)
library(coda)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(cowplot)
library(gtools)

set.seed(123)

# --------------------------
# 1) DATASET
# --------------------------

nb_clusters = 3

load("dataset/soft_connectivity.data.RData")

X = cbind(soft_connectivity.data$longitude_breed,
          soft_connectivity.data$latitude_breed,
          soft_connectivity.data$longitude_wint,
          soft_connectivity.data$latitude_wint) |> 
  as.matrix() 

X_tbl = X |> 
  as.data.frame()

# --------------------------
# 2) CODES NIMBLE
# --------------------------

#  Connectivité seule
code_connectivity_only <- nimbleCode({
  pi[1:K] ~ ddirch(alpha[1:K])
  for(k in 1:K){
    for(j in 1:d){
      mu[k,j] ~ dnorm(0, var=10)
      sigma[k,j] ~ dinvgamma(shape=2.1, scale=1.1)
      prec[k,j] <- 1/sigma[k,j]
    }
  }
  for(i in 1:N){
    cl[i] ~ dcat(pi[1:K])
    for(j in 1:d) X[i,j] ~ dnorm(mu[cl[i],j], prec[cl[i],j])
  }
})

# --------------------------
# 3) CONSTANTES, DATA, INITS
# --------------------------

N = nrow(X)
K=nb_clusters
d=4

d1.mean = mean(X_tbl %>% pull(1), na.rm = TRUE)
d2.mean = mean(X_tbl %>% pull(2), na.rm = TRUE)
d3.mean = mean(X_tbl %>% pull(3), na.rm = TRUE)
d4.mean = mean(X_tbl %>% pull(4), na.rm = TRUE)

mean.position = c(d1.mean, d2.mean, d3.mean, d4.mean)

init.mu1 = mean.position + rnorm(length(mean.position), mean = 0, sd = 1)
init.mu2 = mean.position + rnorm(length(mean.position), mean = 0, sd = 1)
init.mu3 = mean.position + rnorm(length(mean.position), mean = 0, sd = 1)

constants <- list(K=K, N=N, d=d, alpha=rep(1,K), init=c(1,0))
data_conn_only <- list(X = X)


inits_common <- list(
  cl = sample(1:K, N, replace = TRUE),   
  mu = matrix(c(init.mu1, init.mu2, init.mu3), nrow = K, byrow = TRUE),
  sigma = matrix(1, nrow = K, ncol = d)
)


inits_conn <- list(
  cl = inits_common$cl,   
  mu = inits_common$mu,
  sigma = inits_common$sigma,
  pi = rep(1/K, K)
)



# --------------------------
# 4) RUN NIMBLE helper
# --------------------------
run_nimble <- function(code, constants, data, inits, monitors,
                       niter = 10000, nburnin = 2000, thin = 4, progressBar=TRUE) {
  model <- nimbleModel(code, constants = constants, data = data, inits = inits, calculate = TRUE)
  cmodel <- compileNimble(model)
  conf <- configureMCMC(model, monitors = monitors)
  mcmc <- buildMCMC(conf)
  cmcmc <- compileNimble(mcmc, project = model)
  samp <- runMCMC(cmcmc, niter = niter, nburnin = nburnin, thin = thin, progressBar = progressBar)
  return(list(model = model, cmodel = cmodel, samples = as.mcmc(samp)))
}

# --------------------------
# 5) LANCER LES CHAÎNES (ajuster niter si nécessaire)
# --------------------------
niter <- 6000; nburnin <- 2000; thin <- 4


cat("Running connectivity-only model...\n")
res_conn <- run_nimble(code_connectivity_only, constants, data_conn_only, inits_conn,
                       monitors = c("pi","mu","sigma"),
                       niter=niter, nburnin=nburnin, thin=thin)

# --------------------------
# 6) Convertir samples en data.frame
# --------------------------
samples_conn_df <- as.data.frame(as.matrix(res_conn$samples))

samples_conn_df[1:3,]

library(ggplot2)
library(dplyr)
library(purrr)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(cowplot)

# --------------------------------------------------
# 1) Extraire mu_post et Sigma_post (diagonaux)
# --------------------------------------------------
mu_cols <- grep("^mu\\[", colnames(samples_conn_df), value = TRUE)
k_idx <- as.integer(sub("^mu\\[(\\d+),\\s*(\\d+)\\]$", "\\1", mu_cols))
d_idx <- as.integer(sub("^mu\\[(\\d+),\\s*(\\d+)\\]$", "\\2", mu_cols))

K <- max(k_idx); D <- max(d_idx)

mu_post <- matrix(NA_real_, nrow = K, ncol = D)
for (c in seq_along(mu_cols)) {
  mu_post[k_idx[c], d_idx[c]] <- mean(samples_conn_df[[mu_cols[c]]])
}
colnames(mu_post) <- c("lon_breed","lat_breed","lon_wint","lat_wint")

# matrices diagonales
Sigma_post <- vector("list", K)
for (k in 1:K) {
  sigmas <- numeric(D)
  for (d in 1:D) {
    colname <- paste0("sigma[", k, ", ", d, "]")
    sigmas[d] <- mean(samples_conn_df[[colname]])
  }
  Sigma_post[[k]] <- diag(sigmas, D, D)
}

# --------------------------------------------------
# 2) Fonctions ellipses/aplats
# --------------------------------------------------
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

# --------------------------------------------------
# 3) Construire projections
# --------------------------------------------------
cluster_cols <- c("red3", "forestgreen", "royalblue")

# Été = dims 1-2
mu_summer    <- mu_post[, 1:2, drop = FALSE]
Sigma_summer <- lapply(Sigma_post, function(S) S[1:2, 1:2])
polys_summer <- make_aplats(split(mu_summer, row(mu_summer)),
                            Sigma_summer, fill_cols = cluster_cols)
df_cent_summer <- data.frame(x = mu_summer[,1], y = mu_summer[,2], comp = factor(1:K))

# Hiver = dims 3-4
mu_winter    <- mu_post[, 3:4, drop = FALSE]
Sigma_winter <- lapply(Sigma_post, function(S) S[3:4, 3:4])
polys_winter <- make_aplats(split(mu_winter, row(mu_winter)),
                            Sigma_winter, fill_cols = cluster_cols)
df_cent_winter <- data.frame(x = mu_winter[,1], y = mu_winter[,2], comp = factor(1:K))

# --------------------------------------------------
# 4) Carte Europe
# --------------------------------------------------
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
  labs(title = "Barycentres et ellipses - Été",
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
  labs(title = "Barycentres et ellipses - Hiver",
       x = "Longitude", y = "Latitude")

plot_grid(p_summer, p_winter, ncol = 2)

