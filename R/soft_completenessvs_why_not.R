# ==================================================
# GMM Bayésien (NIMBLE) en 4D (été+"hiver") avec correction par Ω
# - Simulation 4D bloc-diagonale: (x_su,y_su) ⟂ (x_wi,y_wi)
# - Les K composantes partagent LE MÊME centre/covariance en hiver (simulation)
#   => l'hiver n'apporte pas d'info de clusterisation
# - Mais le modèle N'IMPOSE PAS cette contrainte (il ne le sait pas)
# - Deux Ω distincts: Ω_su(x,y) et Ω_wi(x,y)
# - LogZ calculés séparément par saison en 2D et utilisés selon la saison d'observation
# - Ones-trick stable et borné comme dans Script 2
# - Modèle naïf en 4D (sans Ω) pour comparaison
# ==================================================

rm(list = ls())
set.seed(123)

suppressPackageStartupMessages({
  library(nimble)
  library(MASS)
  library(Matrix)
  library(terra)
  library(mvtnorm)
  library(ggplot2)
})

# ----------------------- utilitaires -----------------------
make_cov2D <- function(angle_deg, var_major, var_minor){
  phi <- angle_deg * pi/180
  R <- matrix(c(cos(phi), -sin(phi), sin(phi), cos(phi)), 2, 2, byrow = TRUE)
  S <- R %*% diag(c(var_major, var_minor)) %*% t(R)
  ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  if (!isSymmetric(S) || any(ev <= 1e-10)) S <- as.matrix(Matrix::nearPD(S)$mat)
  S
}

gauss2d <- function(x, y, mu, Sigma) {
  invS <- solve(Sigma)
  dx <- x - mu[1]; dy <- y - mu[2]
  Q <- invS[1,1]*dx*dx + (invS[1,2]+invS[2,1])*dx*dy + invS[2,2]*dy*dy
  exp(-0.5 * Q)
}

# ======================= 1) Simulation 4D =======================
K <- 3L; D <- 4L; N_tot <- 6000L
pi_true <- c(0.20, 0.35, 0.45)

# Été (par k)
mu_su_true <- list(
  c(-3.4,  0.0),
  c( 0.0, -0.6),
  c( 3.4,  0.0)
)
Sigma_su_true <- list(
  matrix(c(0.9, 0.0,
           0.0, 3.5), 2, 2, byrow = TRUE),
  matrix(c(6.0, -3.8,
           -3.8,  4.8), 2, 2, byrow = TRUE),  # diagonale ~ -45°
  matrix(c(0.8, 0.0,
           0.0, 3.2), 2, 2, byrow = TRUE)
)

# Hiver (COMMUN à tous les k à la SIMULATION)
mu_wi_true <- c(3.0, -26.0)
Sigma_wi_true <- make_cov2D(angle_deg = 1, var_major = 2.0, var_minor = 2.0)

# Assemble mu/Sigma 4D bloc-diagonaux
mu_true <- vector("list", K)
Sigma_true <- vector("list", K)
for(k in 1:K){
  mu_true[[k]] <- c(mu_su_true[[k]], mu_wi_true)
  Sigma_true[[k]] <- as.matrix(Matrix::bdiag(Sigma_su_true[[k]], Sigma_wi_true))
}

# Echantillonnage latent + données 4D
z <- sample(1:K, N_tot, replace = TRUE, prob = pi_true)
X4_all <- t(vapply(1:N_tot, function(i) MASS::mvrnorm(1, mu = mu_true[[z[i]]], Sigma = Sigma_true[[z[i]]]), numeric(D)))
colnames(X4_all) <- c("x_su","y_su","x_wi","y_wi")
X_su_all <- X4_all[,1:2, drop = FALSE]
X_wi_all <- X4_all[,3:4, drop = FALSE]

# ======================= 2) Ω été / Ω hiver =======================
# Grilles / rasters séparés pour l'intégration 2D
pad <- 1
xmin_all <- min(c(X_su_all[,1], X_wi_all[,1])) - pad
xmax_all <- max(c(X_su_all[,1], X_wi_all[,1])) + pad
ymin_all <- min(c(X_su_all[,2], X_wi_all[,2])) - pad
ymax_all <- max(c(X_su_all[,2], X_wi_all[,2])) + pad

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
  r[] <- with(xy, omega_fun(x, y)) * 0.5
  r
}

build_omega_winter <- function(xmin, xmax, ymin, ymax, nx = 60, ny = 60){
  r <- rast(nrows = ny, ncols = nx, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
  xy <- as.data.frame(xyFromCell(r, 1:ncell(r))); colnames(xy) <- c("x","y")
  center_wi_1 <- c(3, -26)
  Sigma_wi_1  <- make_cov2D(angle_deg =  1, var_major = 2.0, var_minor = 2.0)
  omega_fun <- function(x, y){
    base_wi <- 0.04; cap_wi <- 0.90; amp1 <- 1.00
    val <- base_wi + amp1 * gauss2d(x, y, center_wi_1, Sigma_wi_1)
    pmax(0, pmin(val, cap_wi))
  }
  r[] <- with(xy, omega_fun(x, y)) * 0.3
  r
}

r_su <- build_omega_summer(xmin_all, xmax_all, ymin_all, ymax_all, nx = 60, ny = 60)
r_wi <- build_omega_winter(xmin_all, xmax_all, ymin_all, ymax_all, nx = 60, ny = 60)

Omega_su_all <- terra::extract(r_su, X_su_all)[,1]
Omega_wi_all <- terra::extract(r_wi, X_wi_all)[,1]
Omega_su_all[!is.finite(Omega_su_all)] <- 0
Omega_wi_all[!is.finite(Omega_wi_all)] <- 0

# ======================= 3) Filtrage par Ω et table long =======================
# Chaque individu peut être observé en été OU en hiver (comme Script 1)
sel_ete   <- runif(N_tot) < Omega_su_all
sel_hiver <- runif(N_tot) < Omega_wi_all

obs_ete <- data.frame(
  id = which(sel_ete), saison = "ete",
  x_su = X_su_all[sel_ete,1], y_su = X_su_all[sel_ete,2],
  x_wi = X_wi_all[sel_ete,1], y_wi = X_wi_all[sel_ete,2],
  Omega = Omega_su_all[sel_ete]
)
obs_hiver <- data.frame(
  id = which(sel_hiver), saison = "hiver",
  x_su = X_su_all[sel_hiver,1], y_su = X_su_all[sel_hiver,2],
  x_wi = X_wi_all[sel_hiver,1], y_wi = X_wi_all[sel_hiver,2],
  Omega = Omega_wi_all[sel_hiver]
)

obs_long <- rbind(obs_ete, obs_hiver)
N <- nrow(obs_long)
X <- as.matrix(obs_long[,c("x_su","y_su","x_wi","y_wi")])
ss <- ifelse(obs_long$saison=="ete", 1L, 2L)
Omega_vec <- obs_long$Omega

# ======================= 4) Grilles pour Z_ete, Z_hiver =======================
L <- 60L
ext_all <- ext(r_su)
xmin <- ext_all[1]; xmax <- ext_all[2]
ymin <- ext_all[3]; ymax <- ext_all[4]
dx <- (xmax - xmin) / L; dy <- (ymax - ymin) / L
A_cell <- dx * dy
xs <- xmin + (0:(L-1) + 0.5) * dx
ys <- ymin + (0:(L-1) + 0.5) * dy

grid_mat <- as.matrix(expand.grid(x = xs, y = ys))
omega_su_grid <- as.numeric(terra::extract(r_su, grid_mat, method = "bilinear")[,1])
omega_wi_grid <- as.numeric(terra::extract(r_wi, grid_mat, method = "bilinear")[,1])
omega_su_grid[!is.finite(omega_su_grid)] <- 0
omega_wi_grid[!is.finite(omega_wi_grid)] <- 0
M <- nrow(grid_mat)

# ======================= 5) Fonctions NIMBLE =======================
# MVN (Cholesky)
dmvnorm_nimble <- nimbleFunction(
  run = function(x = double(1), mean = double(1), Prec = double(2),
                 log = logical(0, default = FALSE)) {
    returnType(double(0))
    Dloc <- length(x); xm <- x - mean
    qf <- inprod(xm, Prec %*% xm)
    U <- chol(Prec); ldet <- 2 * sum(log(diag(U)))
    logdens <- 0.5 * ldet - 0.5 * Dloc * log(2*pi) - 0.5 * qf
    if (log) return(logdens) else return(exp(logdens))
  }
)

# log Z saisonnier (2D) par quadrature
logZ_season <- nimbleFunction(
  run=function(mu=double(2), Prec=double(3), pi=double(1),
               grid=double(2), omega=double(1), A=double(0),
               K=integer(0), d1=integer(0), d2=integer(0), M=integer(0)){
    returnType(double(0))
    mean2 <- numeric(2)
    Prec2 <- matrix(0.0, 2, 2)
    x2 <- numeric(2)
    sumZ <- 0.0
    for(m in 1:M){
      x2[1] <- grid[m,1]; x2[2] <- grid[m,2]
      mix <- 0.0
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
    Z <- A * sumZ
    if (Z < 1e-300) Z <- 1e-300
    return(log(Z))
  }
)

# ======================= 6) Modèle corrigé Ω (4D) =======================
C_ub   <- 200
p_min  <- 1e-300
p_max1 <- 1 - 1e-12

code_corr4D <- nimbleCode({
  # Z par saison (2D chacun)
  logZ_su <- logZ_season(mu[1:K,1:D], Prec[1:D,1:D,1:K], pi[1:K],
                         grid[1:M,1:2], omega_su[1:M], A, K, 1, 2, M)
  logZ_wi <- logZ_season(mu[1:K,1:D], Prec[1:D,1:D,1:K], pi[1:K],
                         grid[1:M,1:2], omega_wi[1:M], A, K, 3, 4, M)
  logZ[1] <- logZ_su
  logZ[2] <- logZ_wi
  
  for (i in 1:N) {
    for (k in 1:K) {
      dens[i,k] <- pi[k] * exp(dmvnorm_nimble(X[i,1:D], mu[k,1:D], Prec[1:D,1:D,k], TRUE))
    }
    mixdens[i] <- sum(dens[i,1:K])
    
    ll[i] <- log(Omega[i] * mixdens[i] + 1e-300) - logZ[ss[i]]
    
    p_raw[i]     <- exp(ll[i] - C_ub)
    p_clip_hi[i] <- min(p_raw[i], p_max1)
    p[i]         <- max(p_clip_hi[i], p_min)
    ones[i] ~ dbern(p[i])
  }
  
  # Priors génériques (le modèle ne sait pas que l'hiver est commun)
  pi[1:K] ~ ddirch(alpha[1:K])
  for (k in 1:K) {
    mu[k,1:D] ~ dmnorm(mu0[1:D], prec = Prec0[1:D,1:D])
  }
  
  # Bloc-diagonal (Précisions été/hiver séparées par k, sans contrainte)
  for (k in 1:K) {
    Prec_su[1:2,1:2,k] ~ dwish(R_su[1:2,1:2], df_su)
    Prec_wi[1:2,1:2,k] ~ dwish(R_wi[1:2,1:2], df_wi)
    
    for(i in 1:2) for(j in 1:2) Prec[i,   j,   k] <- Prec_su[i,j,k]
    for(i in 1:2) for(j in 1:2) Prec[i+2, j+2, k] <- Prec_wi[i,j,k]
    for(i in 1:2) for(j in 3:4) Prec[i,   j,   k] <- 0
    for(i in 3:4) for(j in 1:2) Prec[i,   j,   k] <- 0
  }
})

# ======================= 7) Données/constantes =======================
constants_corr <- list(
  N=N, D=D, K=K,
  M=M, A=A_cell,
  alpha=rep(1,K),
  mu0=rep(0,D),
  Prec0=diag(D)*1e-2,
  R_su=diag(2), R_wi=diag(2),
  df_su=4, df_wi=4,
  # ones-trick bounds
  p_min=p_min, p_max1=p_max1, C_ub=C_ub,
  ss=as.integer(ss)
)

data_corr <- list(
  X=X,
  Omega=Omega_vec,
  ones=rep(1L,N),
  grid=grid_mat,
  omega_su=omega_su_grid,
  omega_wi=omega_wi_grid
)

# ======================= 8) Inits (EM facultatif) =======================
# Inits simples et stables
mu_init <- cbind(
  rnorm(K, 0, 0.5) + c(-3, 0, 3),
  rnorm(K, 0, 0.5) + c( 0,-1, 0),
  rnorm(K, 0, 0.2) + mu_wi_true[1],
  rnorm(K, 0, 0.2) + mu_wi_true[2]
)
Prec_su_init <- array(NA_real_, dim = c(2,2,K))
Prec_wi_init <- array(NA_real_, dim = c(2,2,K))
for(k in 1:K){
  Prec_su_init[,,k] <- solve(Sigma_su_true[[k]])
  Prec_wi_init[,,k] <- solve(Sigma_wi_true)
}

inits1 <- list(
  mu = mu_init + matrix(rnorm(K*D, 0, 0.15), K, D),
  pi = rep(1/K, K),
  Prec_su = Prec_su_init,
  Prec_wi = Prec_wi_init
)

# ======================= 9) Build & Run (corrigé Ω) =======================
model4 <- nimbleModel(code_corr4D, data = data_corr, constants = constants_corr, inits = inits1, check = FALSE)
cmodel4 <- compileNimble(model4)
conf4   <- configureMCMC(model4, monitors = c("mu","pi","Prec","logZ"))
mcmc4   <- buildMCMC(conf4)
cmcmc4  <- compileNimble(mcmc4, project = cmodel4)
samples_corr4 <- runMCMC(cmcmc4, niter = 3000, nburnin = 1000, thin = 5, setSeed = TRUE)

cat("\n[Corrigé Ω 4D] Terminé.\n")

# ======================= 10) Modèle naïf 4D (sans Ω) =======================
code_naif4D <- nimbleCode({
  for (i in 1:N){
    z[i] ~ dcat(pi[1:K])
    X[i,1:D] ~ dmnorm(mu[z[i],1:D], prec = Prec[1:D,1:D,z[i]])
  }
  pi[1:K] ~ ddirch(alpha[1:K])
  for (k in 1:K){
    mu[k,1:D] ~ dmnorm(mu0[1:D], prec = Prec0[1:D,1:D])
    Prec[1:D,1:D,k] ~ dwish(R4[1:D,1:D], df4)
  }
})
constants_naif <- list(N=N,D=D,K=K,alpha=rep(1,K),mu0=rep(0,D),Prec0=diag(D)*1e-2,R4=diag(D),df4=D+2)
data_naif <- list(X=X)
inits_naif <- list(z=sample(1:K,N,TRUE), mu = mu_init, Prec = array(0, dim=c(D,D,K)), pi=rep(1/K,K))
for(k in 1:K){ inits_naif$Prec[,,k] <- diag(D) }

model_n4  <- nimbleModel(code_naif4D, data=data_naif, constants=constants_naif, inits=inits_naif, check=FALSE)
cmodel_n4 <- compileNimble(model_n4)
conf_n4   <- configureMCMC(model_n4, monitors=c("mu","pi","Prec"))
mcmc_n4   <- buildMCMC(conf_n4)
cmcmc_n4  <- compileNimble(mcmc_n4, project=model_n4)
samples_naif4 <- runMCMC(cmcmc_n4, niter=4000, nburnin=2000, thin=8, setSeed=TRUE)

cat("[Naïf 4D] Terminé.\n")

# ======================= 11) Petit sanity check =======================
get_pi_means <- function(samples, K){
  pi_cols <- grep("^pi\\[", colnames(samples), value = TRUE)
  ord <- order(as.integer(sub("^pi\\[(\\d+)\\]$", "\\1", pi_cols)))
  colMeans(samples[, pi_cols[ord], drop = FALSE])[1:K]
}

pi_c <- get_pi_means(samples_corr4, K)
pi_n <- get_pi_means(samples_naif4, K)
print(rbind(true=pi_true, corr=round(pi_c,3), naif=round(pi_n,3)))

# Note: pour des comparaisons plus poussées (ellipses, densités),
# on peut reprendre les outils du Script 2 et afficher les marges été/hiver.



