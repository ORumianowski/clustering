# ==================================================
# Modèle NIMBLE propre (GMM 4D + zeros-trick + Ω)
# - pmax() → max() (fix compilation)
# - ss passé dans constants (supprime le warning "non-constant indexes")
# - Grille LxL appariée pour normalisation
# ==================================================

set.seed(123)
suppressPackageStartupMessages({
  library(MASS)
  library(Matrix)   # nearPD
  library(terra)
  library(nimble)
})

# -------------------- utilitaires --------------------
gauss2d <- function(x, y, mu, Sigma) {
  invS <- solve(Sigma)
  dx <- x - mu[1]; dy <- y - mu[2]
  Q <- invS[1,1]*dx*dx + (invS[1,2]+invS[2,1])*dx*dy + invS[2,2]*dy*dy
  exp(-0.5 * Q)
}
build_sigma <- function(sd, R){
  D4 <- diag(sd, 4, 4)
  S <- D4 %*% R %*% D4
  ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  if (!isSymmetric(S) || any(ev <= 1e-10)) S <- as.matrix(nearPD(S)$mat)
  S
}
make_cov2D <- function(angle_deg, var_major, var_minor){
  phi <- angle_deg * pi/180
  R <- matrix(c(cos(phi), -sin(phi), sin(phi), cos(phi)), 2, 2, byrow = TRUE)
  S <- R %*% diag(c(var_major, var_minor)) %*% t(R)
  ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  if (!isSymmetric(S) || any(ev <= 1e-10)) S <- as.matrix(nearPD(S)$mat)
  S
}

# ==================================================
# 1) Paramètres GMM 4D + Simulation
# ==================================================
K <- 3
D <- 4
N_tot <- 3000
pi_true <- c(0.20, 0.45, 0.35)

mu <- list(
  c(-3.5,  0.0,  -3.0,  1.2),
  c( 0.0, -0.6,   0.2, -0.4),
  c( 3.5,  0.0,   3.2,  0.3)
)
for(k in 1:K) mu[[k]][4] <- mu[[k]][4] - 30

sd1 <- c(0.7, 2.8, 2.6, 1.9)
R1 <- matrix(c(
  1.00,  0.00,  0.55, -0.15,
  0.00,  1.00,  0.10,  0.30,
  0.55,  0.10,  1.00,  0.65,
  -0.15,  0.30,  0.65,  1.00
), 4, 4, byrow = TRUE)
Sigma1 <- build_sigma(sd1, R1)

sd2 <- c(2.6, 2.4, 0.8, 3.0)
R2 <- matrix(c(
  1.00, -0.70,  0.40, -0.35,
  -0.70,  1.00, -0.10,  0.45,
  0.40, -0.10,  1.00,  0.00,
  -0.35,  0.45,  0.00,  1.00
), 4, 4, byrow = TRUE)
Sigma2 <- build_sigma(sd2, R2)

sd3 <- c(0.8, 3.2, 3.1, 0.9)
R3 <- matrix(c(
  1.00,  0.00,  0.50,  0.05,
  0.00,  1.00, -0.20,  0.25,
  0.50, -0.20,  1.00,  0.05,
  0.05,  0.25,  0.05,  1.00
), 4, 4, byrow = TRUE)
Sigma3 <- build_sigma(sd3, R3)
Sigma <- list(Sigma1, Sigma2, Sigma3)

z <- sample(1:K, N_tot, replace = TRUE, prob = pi_true)
X4_all <- t(vapply(
  1:N_tot, function(i) MASS::mvrnorm(1, mu = mu[[z[i]]], Sigma = Sigma[[z[i]]]),
  numeric(4)
))
colnames(X4_all) <- c("x_su","y_su","x_wi","y_wi")
X_su_all <- X4_all[, 1:2, drop = FALSE]
X_wi_all <- X4_all[, 3:4, drop = FALSE]

# ==================================================
# 2) Ω(x, saison) + Filtrage
# ==================================================
build_omega_summer <- function(xmin, xmax, ymin, ymax, nx = 60, ny = 60){
  r <- rast(nrows = ny, ncols = nx, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
  xy <- as.data.frame(xyFromCell(r, 1:ncell(r))); colnames(xy) <- c("x","y")
  omega_fun <- function(x, y){
    base <- 0.04
    bumpL <- 0.90 * gauss2d(x, y, c(-3.8, 2.6), diag(c(0.60, 0.60)))
    bumpR <- 0.90 * gauss2d(x, y, c( 3.8, 2.6), diag(c(0.65, 0.65)))
    bumpC <- 1.00 * gauss2d(x, y, c( 0.0,-3.2), diag(c(6.0, 0.9)))
    pmax(0, pmin(base + bumpL + bumpR + bumpC, 0.70))
  }
  r[] <- with(xy, omega_fun(x, y))
  r
}
build_omega_winter <- function(xmin, xmax, ymin, ymax, nx = 60, ny = 60, mu, enlarge = 1.0){
  center_wi_1 <- c(mu[[1]][3], mu[[1]][4]) + c( 0.6, -1.5)
  Sigma_wi_1  <- make_cov2D(angle_deg =  1, var_major = 5.0*enlarge, var_minor = 2.0*enlarge)
  base_wi <- 0.04; cap_wi <- 0.90; amp1 <- 1.00
  r <- rast(nrows = ny, ncols = nx, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
  xy <- as.data.frame(xyFromCell(r, 1:ncell(r))); colnames(xy) <- c("x","y")
  omega_fun <- function(x, y){
    val <- base_wi
    val <- val + amp1 * gauss2d(x, y, center_wi_1, Sigma_wi_1)
    pmax(0, pmin(val, cap_wi))
  }
  r[] <- with(xy, omega_fun(x, y))
  r
}
pad <- 1
xmin_all <- min(c(X_su_all[,1], X_wi_all[,1])) - pad
xmax_all <- max(c(X_su_all[,1], X_wi_all[,1])) + pad
ymin_all <- min(c(X_su_all[,2], X_wi_all[,2])) - pad
ymax_all <- max(c(X_su_all[,2], X_wi_all[,2])) + pad
r_su <- build_omega_summer(xmin_all, xmax_all, ymin_all, ymax_all, nx = 60, ny = 60)
r_wi <- build_omega_winter(xmin_all, xmax_all, ymin_all, ymax_all, nx = 60, ny = 60, mu = mu)

Omega_su_all <- terra::extract(r_su, X_su_all)[,1]
Omega_wi_all <- terra::extract(r_wi, X_wi_all)[,1]

sel_ete  <- runif(N_tot) < Omega_su_all
sel_hiver<- runif(N_tot) < Omega_wi_all
obs_ete  <- data.frame(id=which(sel_ete), saison="été",
                       x_su=X_su_all[sel_ete,1], y_su=X_su_all[sel_ete,2],
                       x_wi=X_wi_all[sel_ete,1], y_wi=X_wi_all[sel_ete,2],
                       Omega=Omega_su_all[sel_ete])
obs_hiver<- data.frame(id=which(sel_hiver), saison="hiver",
                       x_su=X_su_all[sel_hiver,1], y_su=X_su_all[sel_hiver,2],
                       x_wi=X_wi_all[sel_hiver,1], y_wi=X_wi_all[sel_hiver,2],
                       Omega=Omega_wi_all[sel_hiver])
obs_long <- rbind(obs_ete, obs_hiver)
N <- nrow(obs_long)

# ==================================================
# 3) Grille LxL appariée (intégration)
# ==================================================
D <- 4L
X <- as.matrix(obs_long[,c("x_su","y_su","x_wi","y_wi")])
Omega <- obs_long$Omega
ss <- ifelse(obs_long$saison=="été",1L,2L)  # -> sera mis dans constants
zeros <- rep(0L,N)

L <- 60L
ext_all <- ext(r_su)  # même emprise que r_wi
xmin <- ext_all[1]; xmax <- ext_all[2]
ymin <- ext_all[3]; ymax <- ext_all[4]
dx <- (xmax - xmin) / L
dy <- (ymax - ymin) / L
A_cell <- dx * dy
xs <- xmin + (0:(L-1) + 0.5) * dx
ys <- ymin + (0:(L-1) + 0.5) * dy
grid_mat <- as.matrix(expand.grid(x = xs, y = ys))  # M = L^2

omega_su_grid <- as.numeric(terra::extract(r_su, grid_mat, method = "bilinear")[,1])
omega_wi_grid <- as.numeric(terra::extract(r_wi, grid_mat, method = "bilinear")[,1])

grid_su <- grid_mat
grid_wi <- grid_mat
M_su <- nrow(grid_su)
M_wi <- nrow(grid_wi)
A_su <- A_cell
A_wi <- A_cell

# ==================================================
# 4) NIMBLE : fonctions + modèle
# ==================================================
dmvnorm_nimble <- nimbleFunction(
  run=function(x=double(1), mean=double(1), Prec=double(2),
               log=logical(0, default=FALSE)){
    returnType(double(0))
    xm <- x - mean
    qf <- inprod(xm, Prec %*% xm)
    U <- chol(Prec)
    ldet <- 2 * sum(log(diag(U)))
    logdens <- 0.5*ldet - 0.5*length(x)*log(2*pi) - 0.5*qf
    if(log) return(logdens) else return(exp(logdens))
  }
)

logZ_calc <- nimbleFunction(
  run=function(mu=double(2), Prec=double(3), pi=double(1),
               grid=double(2), omega=double(1), A=double(0),
               K=integer(0), Dfull=integer(0),
               d1=integer(0), d2=integer(0), M=integer(0)){
    returnType(double(0))
    mean2 <- numeric(2)
    Prec2 <- matrix(0.0, 2, 2)
    x2 <- numeric(2)
    sumZ <- 0.0
    for(m in 1:M){
      x2[1] <- grid[m,1]; x2[2] <- grid[m,2]
      mix <- 0
      for(k in 1:K){
        mean2[1] <- mu[k,d1]; mean2[2] <- mu[k,d2]
        Prec2[1,1] <- Prec[d1,d1,k]
        Prec2[1,2] <- Prec[d1,d2,k]
        Prec2[2,1] <- Prec[d2,d1,k]
        Prec2[2,2] <- Prec[d2,d2,k]
        mix <- mix + pi[k]*dmvnorm_nimble(x2, mean2, Prec2, FALSE)
      }
      sumZ <- sumZ + omega[m]*mix
    }
    Z <- A*(sumZ/M)
    if(Z < 1e-300) Z <- 1e-300
    return(log(Z))
  }
)

code_corrige <- nimbleCode({
  logZ_su <- logZ_calc(mu[1:K,1:D], Prec[1:D,1:D,1:K], pi[1:K],
                       grid_su[1:M_su,1:2], omega_su_grid[1:M_su],
                       A_su, K, D, 1, 2, M_su)
  logZ_wi <- logZ_calc(mu[1:K,1:D], Prec[1:D,1:D,1:K], pi[1:K],
                       grid_wi[1:M_wi,1:2], omega_wi_grid[1:M_wi],
                       A_wi, K, D, 3, 4, M_wi)
  logZ[1] <- logZ_su
  logZ[2] <- logZ_wi
  
  for(i in 1:N){
    for(k in 1:K){
      dens[i,k] <- pi[k] * exp(dmvnorm_nimble(X[i,1:D], mu[k,1:D],
                                              Prec[1:D,1:D,k], TRUE))
    }
    mixdens[i] <- sum(dens[i,1:K])
    ll[i] <- log(Omega[i]*mixdens[i] + 1e-300) - logZ[ss[i]]
    ll_clip[i] <- min(max(ll[i], -700), 700)
    zeros[i] ~ dpois(phi[i])
    phi[i] <- max(-ll_clip[i] + 1e-6, 1e-6)  # <- FIX pmax → max
  }
  
  pi[1:K] ~ ddirch(alpha[1:K])
  for(k in 1:K){
    mu[k,1:D] ~ dmnorm(mu0[1:D], prec=Prec0[1:D,1:D])
    Prec[1:D,1:D,k] ~ dwish(R[1:D,1:D], df)
  }
})

# ==================================================
# 5) Données/constantes pour NIMBLE
#    (ss passe dans constants → plus de warning "non-constant indexes")
# ==================================================
constants <- list(
  N=N, D=D, K=K,
  M_su=M_su, M_wi=M_wi,
  A_su=A_su, A_wi=A_wi,
  alpha=rep(1,K), mu0=rep(0,D),
  Prec0=diag(D)*1e-2, R=diag(D), df=D+2,
  ss=as.integer(ss)  # <- ICI dans constants
)

data_list <- list(
  X=X,
  Omega=Omega,
  zeros=as.integer(zeros),
  grid_su=grid_su, grid_wi=grid_wi,
  omega_su_grid=omega_su_grid, omega_wi_grid=omega_wi_grid
)

# ==================================================
# 6) Inits
# ==================================================
barX <- colMeans(X)
mu_init <- matrix(NA,K,D); for(k in 1:K) mu_init[k,] <- barX + rnorm(D,0,1)
Prec_init <- array(0,dim=c(D,D,K)); for(k in 1:K) Prec_init[,,k] <- diag(D)*0.5
inits1 <- list(mu=mu_init, Prec=Prec_init, pi=rep(1/K,K))
inits2 <- list(mu=mu_init+rnorm(K*D,0,0.5), Prec=Prec_init, pi=rep(1/K,K))

# ==================================================
# 7) Modèle + compile + MCMC (2 chaînes)
# ==================================================
model <- nimbleModel(code_corrige, data = data_list, constants = constants,
                     inits = inits1, check = FALSE)

conf  <- configureMCMC(model, monitors = c("mu","pi","Prec"))
mcmc  <- buildMCMC(conf)

# (1) Compiler le modèle
Cmodel <- nimble::compileNimble(model, showCompilerOutput = TRUE)

# (2) Compiler le MCMC en le reliant au modèle compilé
cmcmc  <- nimble::compileNimble(mcmc, project = Cmodel, showCompilerOutput = TRUE)

# (3) Chaîne 1
cat("➡️ Chaîne 1...\n")
samples1 <- runMCMC(cmcmc, niter = 2500, nburnin = 1000, thin = 2,
                    inits = inits1, setSeed = 100)

# (4) Chaîne 2 (optionnel)
cat("➡️ Chaîne 2...\n")
samples2 <- runMCMC(cmcmc, niter = 2500, nburnin = 1000, thin = 2,
                    inits = inits2, setSeed = 200)

# =============================
# Plots propres pour μ et π
# =============================
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

dim_labels <- c("x_su","y_su","x_wi","y_wi")  # D=4

S1 <- as.data.frame(samples1)
S2 <- as.data.frame(samples2)
S1$.iter <- seq_len(nrow(S1))
S2$.iter <- seq_len(nrow(S2))
S1$.chain <- "chaine 1"
S2$.chain <- "chaine 2"

# ---- helper: jolies étiquettes mu[k,d] et pi[k]
pretty_mu_name <- function(k, d) sprintf("mu[k=%d, %s]", k, dim_labels[d])
pretty_pi_name <- function(k)    sprintf("pi[k=%d]", k)

# ---- construire un long df pour mu
mu_cols <- grep("^mu\\[\\s*\\d+\\s*,\\s*\\d+\\s*\\]$", colnames(S1), value = TRUE)

parse_mu <- function(nm) {
  m <- str_match(nm, "^mu\\[\\s*(\\d+)\\s*,\\s*(\\d+)\\s*\\]$")
  k <- as.integer(m[,2]); d <- as.integer(m[,3])
  tibble(param = nm, k = k, d = d, label = pretty_mu_name(k, d))
}

mu_map <- bind_rows(lapply(mu_cols, parse_mu))

long_mu <- bind_rows(
  S1 %>% select(all_of(mu_cols), .iter, .chain) %>% pivot_longer(cols = all_of(mu_cols), names_to = "param", values_to = "value"),
  S2 %>% select(all_of(mu_cols), .iter, .chain) %>% pivot_longer(cols = all_of(mu_cols), names_to = "param", values_to = "value")
) %>% left_join(mu_map, by = "param")

# ---- construire un long df pour pi
pi_cols <- grep("^pi\\[\\s*\\d+\\s*\\]$", colnames(S1), value = TRUE)

parse_pi <- function(nm) {
  m <- str_match(nm, "^pi\\[\\s*(\\d+)\\s*\\]$")
  k <- as.integer(m[,2])
  tibble(param = nm, k = k, label = pretty_pi_name(k))
}

pi_map <- bind_rows(lapply(pi_cols, parse_pi))

long_pi <- bind_rows(
  S1 %>% select(all_of(pi_cols), .iter, .chain) %>% pivot_longer(cols = all_of(pi_cols), names_to = "param", values_to = "value"),
  S2 %>% select(all_of(pi_cols), .iter, .chain) %>% pivot_longer(cols = all_of(pi_cols), names_to = "param", values_to = "value")
) %>% left_join(pi_map, by = "param")

# ---- thèmes
theme_clean <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        strip.text = element_text(face = "bold"))

# ---- TRACEPLOTS μ
p_trace_mu <- ggplot(long_mu, aes(x = .iter, y = value, color = .chain)) +
  geom_line(alpha = 0.8, linewidth = 0.5) +
  facet_wrap(~ label, scales = "free_y", ncol = 3) +
  labs(title = "Traceplots — μ par composante et dimension",
       x = "itération", y = "valeur", color = "chaîne") +
  theme_clean

# ---- DENSITÉS μ
p_dens_mu <- ggplot(long_mu, aes(x = value, color = .chain, fill = .chain)) +
  geom_density(alpha = 0.25, linewidth = 0.7) +
  facet_wrap(~ label, scales = "free", ncol = 3) +
  labs(title = "Densités postérieures — μ",
       x = "valeur", y = "densité", color = "chaîne", fill = "chaîne") +
  theme_clean

# ---- TRACEPLOTS π
p_trace_pi <- ggplot(long_pi, aes(x = .iter, y = value, color = .chain)) +
  geom_line(alpha = 0.8, linewidth = 0.6) +
  facet_wrap(~ label, scales = "free_y") +
  labs(title = "Traceplots — π (poids des composantes)",
       x = "itération", y = "valeur", color = "chaîne") +
  theme_clean

# ---- DENSITÉS π
p_dens_pi <- ggplot(long_pi, aes(x = value, color = .chain, fill = .chain)) +
  geom_density(alpha = 0.25, linewidth = 0.8) +
  facet_wrap(~ label, scales = "free") +
  labs(title = "Densités postérieures — π",
       x = "valeur", y = "densité", color = "chaîne", fill = "chaîne") +
  theme_clean

# ---- afficher (tu peux les imprimer un par un)
print(p_trace_mu)
print(p_dens_mu)
print(p_trace_pi)
print(p_dens_pi)


# -------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(viridisLite)
  library(mvtnorm)
  library(patchwork)
  library(stringr)
  library(dplyr)
})

# ---------- Fusion des chaînes ----------
samples_corr <- as.data.frame(rbind(samples1, samples2))
have_naif <- exists("samples_naif")

# ---------- Helpers plots ----------
extr_muSigma4D <- function(samples, K, D = 4){
  ## -- mu : parse propre --
  mu_cols <- grep("^mu\\[\\s*\\d+\\s*,\\s*\\d+\\s*\\]$", colnames(samples), value = TRUE)
  k_idx <- as.integer(sub("^mu\\[\\s*(\\d+)\\s*,\\s*(\\d+)\\s*\\]$", "\\1", mu_cols))
  d_idx <- as.integer(sub("^mu\\[\\s*(\\d+)\\s*,\\s*(\\d+)\\s*\\]$", "\\2", mu_cols))
  mu_post <- matrix(NA_real_, nrow = K, ncol = D)
  for (c in seq_along(mu_cols)) mu_post[k_idx[c], d_idx[c]] <- mean(samples[[mu_cols[c]]])
  
  ## -- Prec : parse (i,j,k) et reconstruction sûre --
  Sigma_post <- vector("list", K)
  for (k in 1:K) {
    pat <- sprintf("^Prec\\[\\s*(\\d+)\\s*,\\s*(\\d+)\\s*,\\s*%d\\s*\\]$", k)
    prec_cols <- grep(pat, colnames(samples), value = TRUE)
    if (length(prec_cols) != D*D)
      stop("Colonnes Prec manquantes pour k=", k)
    # parse i,j
    ij <- str_match(prec_cols, "^Prec\\[\\s*(\\d+)\\s*,\\s*(\\d+)\\s*,")[,2:3]
    ii <- as.integer(ij[,1]); jj <- as.integer(ij[,2])
    Prec_k <- matrix(NA_real_, D, D)
    vals <- colMeans(samples[, prec_cols, drop = FALSE])
    for (m in seq_along(prec_cols)) Prec_k[ii[m], jj[m]] <- vals[m]
    # symétrisation douce
    S <- tryCatch(solve(Prec_k), error = function(e) NA)
    if (any(is.na(S))) {
      # petite régularisation si nécessaire
      Prec_k <- Prec_k + diag(D) * 1e-6
      S <- solve(Prec_k)
    }
    S <- (S + t(S)) / 2
    Sigma_post[[k]] <- S
  }
  list(mu = mu_post, Sigma = Sigma_post)
}

ellipse_points <- function(mu, Sigma, r, n = 200) {
  ang <- seq(0, 2*pi, length.out = n)
  circle <- cbind(cos(ang), sin(ang))
  R <- chol(Sigma)
  pts <- circle %*% t(R) * r
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
  pi_cols <- grep("^pi\\[\\s*\\d+\\s*\\]$", colnames(samples), value = TRUE)
  ord <- order(as.integer(sub("^pi\\[\\s*(\\d+)\\s*\\]$", "\\1", pi_cols)))
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
cols_season <- c("été" = "orange", "hiver" = "darkblue")

# ==================================================
# (A) FIGURE points simulés & Ω
# ==================================================
# AVANT Ω (population complète)
before_df <- rbind(
  data.frame(x = X_su_all[,1], y = X_su_all[,2], saison = "été"),
  data.frame(x = X_wi_all[,1], y = X_wi_all[,2], saison = "hiver")
)
p_before <- ggplot(before_df, aes(x = x, y = y, color = saison)) +
  geom_point(size = 0.55, alpha = 0.85) +
  scale_color_manual(values = cols_season) +
  coord_equal(xlim = xlim_all, ylim = ylim_all, expand = FALSE) +
  square_theme + labs(title = "AVANT Ω — Été & Hiver", x = "x", y = "y", color = "Saison")

# Ω ÉTÉ / Ω HIVER
su_r_df <- as.data.frame(r_su, xy = TRUE); names(su_r_df)[3] <- "omega"
wi_r_df <- as.data.frame(r_wi, xy = TRUE); names(wi_r_df)[3] <- "omega"
p_su_omega <- ggplot(su_r_df, aes(x = x, y = y, fill = omega)) +
  geom_raster() +
  coord_equal(xlim = xlim_all, ylim = ylim_all, expand = FALSE) +
  scale_fill_viridis_c(name = "Ω (été)") +
  square_theme + labs(title = "Fond Ω — ÉTÉ", x = "x", y = "y")
p_wi_omega <- ggplot(wi_r_df, aes(x = x, y = y, fill = omega)) +
  geom_raster() +
  coord_equal(xlim = xlim_all, ylim = ylim_all, expand = FALSE) +
  scale_fill_viridis_c(name = "Ω (hiver)") +
  square_theme + labs(title = "Fond Ω — HIVER", x = "x", y = "y")

# APRÈS Ω (observations effectives)
after_df <- rbind(
  transform(obs_long[obs_long$saison == "été", ],
            x = x_su, y = y_su)[, c("x","y","saison")],
  transform(obs_long[obs_long$saison == "hiver", ],
            x = x_wi, y = y_wi)[, c("x","y","saison")]
)
p_after <- ggplot(after_df, aes(x = x, y = y, color = saison)) +
  geom_point(size = 0.55, alpha = 0.9) +
  scale_color_manual(values = cols_season) +
  coord_equal(xlim = xlim_all, ylim = ylim_all, expand = FALSE) +
  square_theme + labs(title = "APRÈS Ω — Été & Hiver", x = "x", y = "y", color = "Saison")

# Disposition (figure "points simulés")
fig_points <- (p_before | p_after) / (p_su_omega | p_wi_omega)

# ==================================================
# (B) Densités séparées Été/Hiver : Vrai vs (Naïf*) vs Corrigé
# ==================================================
K <- length(Sigma); true_Sigma <- Sigma
# Vrai (paramètres de simulation)
mu_true_mat <- do.call(rbind, mu)
mu_su_true  <- mu_true_mat[, 1:2, drop=FALSE]
mu_wi_true  <- mu_true_mat[, 3:4, drop=FALSE]
Sig_su_true <- lapply(true_Sigma, \(S) S[1:2, 1:2, drop = FALSE])
Sig_wi_true <- lapply(true_Sigma, \(S) S[3:4, 3:4, drop = FALSE])
pi_true_vec <- if (exists("pi_true")) pi_true else rep(1/K, K)

# Posterior corrigé
post_corr <- extr_muSigma4D(samples_corr, K, D = 4)
mu_su_corr  <- post_corr$mu[, 1:2, drop=FALSE]
mu_wi_corr  <- post_corr$mu[, 3:4, drop=FALSE]
Sig_su_corr <- lapply(post_corr$Sigma, \(S) S[1:2, 1:2, drop=FALSE])
Sig_wi_corr <- lapply(post_corr$Sigma, \(S) S[3:4, 3:4, drop=FALSE])
pi_corr <- get_pi_means(samples_corr, K)

# (Optionnel) Posterior naïf si disponible
if (have_naif) {
  post_naif <- extr_muSigma4D(samples_naif, K, D = 4)
  mu_su_naif  <- post_naif$mu[, 1:2, drop=FALSE]
  mu_wi_naif  <- post_naif$mu[, 3:4, drop=FALSE]
  Sig_su_naif <- lapply(post_naif$Sigma, \(S) S[1:2, 1:2, drop=FALSE])
  Sig_wi_naif <- lapply(post_naif$Sigma, \(S) S[3:4, 3:4, drop=FALSE])
  pi_naif <- get_pi_means(samples_naif, K)
}

# Grilles pour densités (depuis les rasters)
su_xy <- as.data.frame(r_su, xy = TRUE)[, c("x","y")]
wi_xy <- as.data.frame(r_wi, xy = TRUE)[, c("x","y")]

dens_su_true <- as.data.frame(mix_density_df(mu_su_true, Sig_su_true, pi_true_vec, su_xy))
dens_wi_true <- as.data.frame(mix_density_df(mu_wi_true, Sig_wi_true, pi_true_vec, wi_xy))

dens_su_corr <- as.data.frame(mix_density_df(mu_su_corr, Sig_su_corr, pi_corr, su_xy))
dens_wi_corr <- as.data.frame(mix_density_df(mu_wi_corr, Sig_wi_corr, pi_corr, wi_xy))

if (have_naif) {
  dens_su_naif <- as.data.frame(mix_density_df(mu_su_naif, Sig_su_naif, pi_naif, su_xy))
  dens_wi_naif <- as.data.frame(mix_density_df(mu_wi_naif, Sig_wi_naif, pi_naif, wi_xy))
}

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
df_wi_true_cent <- data.frame(x = mu_wi_true[,1], y = mu_wi_true[,2])
df_su_corr_cent <- data.frame(x = mu_su_corr[,1], y = mu_su_corr[,2])
df_wi_corr_cent <- data.frame(x = mu_wi_corr[,1], y = mu_wi_corr[,2])

p_su_d_true <- p_gmm(dens_su_true, df_su_true_cent, "ÉTÉ — Densité GMM (vrai)",
                     paste0("π = [", paste(sprintf("%.2f", pi_true_vec), collapse = ", "), "]"))
p_wi_d_true <- p_gmm(dens_wi_true, df_wi_true_cent, "HIVER — Densité GMM (vrai)",
                     paste0("π = [", paste(sprintf("%.2f", pi_true_vec), collapse = ", "), "]"))

p_su_d_corr <- p_gmm(dens_su_corr, df_su_corr_cent, "ÉTÉ — Densité GMM (corrigé)",
                     paste0("π = [", paste(sprintf("%.2f", pi_corr), collapse = ", "), "]"))
p_wi_d_corr <- p_gmm(dens_wi_corr, df_wi_corr_cent, "HIVER — Densité GMM (corrigé)",
                     paste0("π = [", paste(sprintf("%.2f", pi_corr), collapse = ", "), "]"))

if (have_naif) {
  df_su_naif_cent <- data.frame(x = mu_su_naif[,1], y = mu_su_naif[,2])
  df_wi_naif_cent <- data.frame(x = mu_wi_naif[,1], y = mu_wi_naif[,2])
  p_su_d_naif <- p_gmm(dens_su_naif, df_su_naif_cent, "ÉTÉ — Densité GMM (naïf)",
                       paste0("π = [", paste(sprintf("%.2f", pi_naif), collapse = ", "), "]"))
  p_wi_d_naif <- p_gmm(dens_wi_naif, df_wi_naif_cent, "HIVER — Densité GMM (naïf)",
                       paste0("π = [", paste(sprintf("%.2f", pi_naif), collapse = ", "), "]"))
}

if (have_naif) {
  row_su <- (p_su_d_true | p_su_d_naif | p_su_d_corr)
  row_wi <- (p_wi_d_true | p_wi_d_naif | p_wi_d_corr)
} else {
  row_su <- (p_su_d_true | p_su_d_corr)
  row_wi <- (p_wi_d_true | p_wi_d_corr)
}
p_top  <- (row_su / row_wi) +
  patchwork::plot_annotation(title = "Densités séparées — Été (haut) / Hiver (bas)") &
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# ==================================================
# (C) Ellipses combinées (été+hiver) + segments par individu
# ==================================================
# sous-ensemble : individus observés au moins une fois
keep_any <- sel_ete | sel_hiver
X_obs_any <- X4_all[keep_any,,drop=FALSE]
z_true_obs <- factor(z[keep_any], levels = as.character(1:K))

# Assignations MAP (corrigé) en 4D
assign_corr <- assign_clusters_MAP(X_obs_any, post_corr$mu, post_corr$Sigma, pi_corr)
z_corr_obs  <- factor(assign_corr, levels = as.character(1:K))

# DataFrames segments été->hiver
pairs_true <- data.frame(lonA = X_obs_any[,1], latA = X_obs_any[,2],
                         lonB = X_obs_any[,3], latB = X_obs_any[,4],
                         cluster = z_true_obs)
pairs_corr <- transform(pairs_true, cluster = z_corr_obs)

# Ellipses combinées (été + hiver) pour Vrai / Corrigé
mk_polys_sets <- function(mu2d_su, Sig2d_su, mu2d_wi, Sig2d_wi){
  make_aplats_sets(
    mu_list_su   = split(mu2d_su, row(mu2d_su)),
    Sigma_list_su= Sig2d_su,
    mu_list_wi   = split(mu2d_wi, row(mu2d_wi)),
    Sigma_list_wi= Sig2d_wi,
    probs = c(0.50, 0.75, 0.90)
  )
}
polys_true <- mk_polys_sets(mu_su_true, Sig_su_true, mu_wi_true, Sig_wi_true)
polys_corr <- mk_polys_sets(mu_su_corr, Sig_su_corr, mu_wi_corr, Sig_wi_corr)

cluster_palette <- setNames(rocket(K, begin = 0.10, end = 0.75, direction = -1),
                            as.character(1:K))
coverage_alphas <- c("50%" = 0.50, "75%" = 0.20, "90%" = 0.10)

build_pairs_plot <- function(pairs_df, polys_df, title){
  ggplot() +
    geom_segment(data = pairs_df,
                 aes(x = lonA, y = latA, xend = lonB, yend = latB, color = cluster),
                 alpha = 0.05, linewidth = 0.5) +
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
p_pairs_corr <- build_pairs_plot(pairs_corr, polys_corr, "Ellipses combinées — CORRIGÉ (été + hiver)")

row_pairs <- (p_pairs_true | p_pairs_corr) +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# ==================================================
# Impression des trois figures
# ==================================================
print(fig_points)
print(p_top)
print(row_pairs)

