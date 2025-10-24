# ==================================================
# GMM bayésien (NIMBLE) en 4D (été + hiver) avec correction par Ω
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

# ----------------------- Utilitaires -----------------------
make_cov2D <- function(angle_deg, var_major, var_minor){
  phi <- angle_deg * pi/180
  R <- matrix(c(cos(phi), -sin(phi), sin(phi), cos(phi)), 2, 2, byrow = TRUE)
  S <- R %*% diag(c(var_major, var_minor)) %*% t(R)
  ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  if (!isSymmetric(S) || any(ev <= 1e-10)) S <- as.matrix(Matrix::nearPD(S)$mat)
  S
}

build_sigma <- function(sd, R){
  D4 <- diag(sd, 4, 4)
  S <- D4 %*% R %*% D4
  ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  if (!isSymmetric(S) || any(ev <= 1e-10)) S <- as.matrix(nearPD(S)$mat)
  S
}

gauss2d <- function(x, y, mu, Sigma) {
  invS <- solve(Sigma)
  dx <- x - mu[1]; dy <- y - mu[2]
  Q <- invS[1,1]*dx*dx + (invS[1,2]+invS[2,1])*dx*dy + invS[2,2]*dy*dy
  exp(-0.5 * Q)
}

# ======================= 1) Simulation 4D =======================
# Nombre de clusters
K <- 3
# Dimension
D <- 4
# Taille de la population
N_tot <- 5000
# Répartition de la population entre les clusters
pi_true <- c(0.20, 0.35, 0.45)

# Centres des clusters
# (x_su, y_su, x_wi, y_wi)
mu_true <- list(
  c(-3.5,  0.0,  -1,   2), 
  c( 0.0, -0.6,  -1,   2),  
  c( 3.5,  0.0,  2,  5)  
)
# Décalage global de la coordonnée y_hiver
for (k in 1:K) mu_true[[k]][4] <- mu_true[[k]][4] - 30

# Bloc-diagonal : (1,2) ⟂ (3,4) — covariances été/hiver nulles
sd1 <- c(0.7, 2.8, 1.0, 2.0)
R1 <- matrix(c(
  1.00,  0.00,  0.00,  0.00,
  0.00,  1.00,  0.00,  0.00,
  0.00,  0.00,  1.00,  0.20,
  0.00,  0.00,  0.20,  1.00
), 4, 4, byrow = TRUE)
Sigma1 <- build_sigma(sd1, R1)

sd2 <- c(2.6, 2.4, 1.0, 2.0)
R2 <- matrix(c(
  1.00, -0.70,  0.00,  0.00,
  -0.70,  1.00,  0.00,  0.00,
  0.00,  0.00,  1.00,  0.20,
  0.00,  0.00,  0.20,  1.00
), 4, 4, byrow = TRUE)
Sigma2 <- build_sigma(sd2, R2)

sd3 <- c(0.8, 3.2, 0.6, 2.8)
R3 <- matrix(c(
  1.00,  0.00,  0.00,  0.00,
  0.00,  1.00,  0.00,  0.00,
  0.00,  0.00,  1.00,  0.00,
  0.00,  0.00,  0.00,  1.00
), 4, 4, byrow = TRUE)
Sigma3 <- build_sigma(sd3, R3)

Sigma_true <- list(Sigma1, Sigma2, Sigma3)

# Échantillonnage latent + données 4D
# Attribution du cluster selon pi
z <- sample(1:K, N_tot, replace = TRUE, prob = pi_true)

# Échantillonnage des coordonnées selon la gaussienne du cluster
X4_all <- t(vapply(1:N_tot, function(i) MASS::mvrnorm(1, mu = mu_true[[z[i]]], Sigma = Sigma_true[[z[i]]]), numeric(D)))

colnames(X4_all) <- c("x_su","y_su","x_wi","y_wi")
X_su_all <- X4_all[,1:2, drop = FALSE]
X_wi_all <- X4_all[,3:4, drop = FALSE]

# ======================= 2) Ω été / Ω hiver =======================
# Grilles/rasters séparés pour l’intégration 2D
pad <- 1
xmin_all <- min(c(X_su_all[,1], X_wi_all[,1])) - pad
xmax_all <- max(c(X_su_all[,1], X_wi_all[,1])) + pad
ymin_all <- min(c(X_su_all[,2], X_wi_all[,2])) - pad
ymax_all <- max(c(X_su_all[,2], X_wi_all[,2])) + pad

# Fonctions de construction de l’effort d’échantillonnage Ω (été et hiver)
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
  center_wi_1 <- c(2, -25)
  Sigma_wi_1  <- make_cov2D(angle_deg =  1, var_major = 2.0, var_minor = 2.0)
  omega_fun <- function(x, y){
    base_wi <- 0.04; cap_wi <- 0.90; amp1 <- 1.00
    val <- base_wi + amp1 * gauss2d(x, y, center_wi_1, Sigma_wi_1)
    pmax(0, pmin(val, cap_wi))
  }
  r[] <- with(xy, omega_fun(x, y)) * 0.4
  r
}

# Nombre de pixels par côté (grille spatiale pour la simulation)
n_pix_side <- 60

# Grilles d’effort d’échantillonnage
r_su <- build_omega_summer(xmin_all, xmax_all, ymin_all, ymax_all, nx = n_pix_side, ny = n_pix_side)
r_wi <- build_omega_winter(xmin_all, xmax_all, ymin_all, ymax_all, nx = n_pix_side, ny = n_pix_side)

# Effort d’échantillonnage associé à chaque individu
Omega_su_all <- terra::extract(r_su, X_su_all)[,1]
Omega_wi_all <- terra::extract(r_wi, X_wi_all)[,1]


# ======================= 3) Filtrage par Ω et table « long » =======================
# Chaque individu peut être observé en été ou en hiver
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
sampling.season <- ifelse(obs_long$saison=="ete", 1L, 2L)
Omega_vec <- obs_long$Omega

# ======================= 4) Grilles pour Z_ete et Z_hiver =======================
# Nombre de pixels par côté pour l’intégration numérique
L <- 60L

# Limites de la zone à intégrer
ext_all <- ext(r_su)
xmin <- ext_all[1]; xmax <- ext_all[2]
ymin <- ext_all[3]; ymax <- ext_all[4]
dx <- (xmax - xmin) / L; dy <- (ymax - ymin) / L

# Aire d’un pixel
A_cell <- dx * dy

# Grille systématique pour l’intégration numérique
xs <- xmin + (0:(L-1) + 0.5) * dx
ys <- ymin + (0:(L-1) + 0.5) * dy
grid_mat <- as.matrix(expand.grid(x = xs, y = ys))

# Valeurs d’effort d’échantillonnage sur la grille
omega_su_grid <- as.numeric(terra::extract(r_su, grid_mat, method = "bilinear")[,1])
omega_wi_grid <- as.numeric(terra::extract(r_wi, grid_mat, method = "bilinear")[,1])

# Nombre total de points de grille
M <- nrow(grid_mat)

# ======================= 5) Fonctions NIMBLE =======================

# dmvnorm_nimble : densité (ou log-densité) d’une normale multivariée
# Écriture en précision (Q = Σ^{-1}) + calculs stables via Cholesky
# Formule (écriture précision) :
#   log f(x) = 0.5*log|Q| - 0.5*d*log(2π) - 0.5*(x-μ)^T Q (x-μ)
dmvnorm_nimble <- nimbleFunction(
  run = function(
    x    = double(1),                 # vecteur d’observation (dimension d)
    mean = double(1),                 # moyenne μ (longueur d)
    Prec = double(2),                 # précision Q = Σ^{-1} (d x d)
    log  = logical(0, default = FALSE)
  ) {
    # Type de retour : scalaire double
    returnType(double(0))
    
    # 1) Pré-calculs
    Dloc <- length(x)                 # d = dimension du vecteur
    xm   <- x - mean                  # centrage : (x - μ)
    
    # 2) Terme quadratique (x-μ)^T Q (x-μ)
    qf <- inprod(xm, Prec %*% xm)
    
    # 3) log|Q| via Cholesky
    U    <- chol(Prec)                # décomposition de Cholesky de Q
    ldet <- 2 * sum(log(diag(U)))     # log|Q|
    
    # 4) Log-densité (écriture précision)
    logdens <- 0.5 * ldet - 0.5 * Dloc * log(2*pi) - 0.5 * qf
    
    if (log) {
      return(logdens)                 # renvoie log f(x)
    } else {
      return(exp(logdens))            # renvoie f(x)
    }
  }
)

# log Z saisonnier (2D) par quadrature sur une grille de M points
logZ_season <- nimbleFunction(
  run=function(mu=double(2),      # K x D : moyennes des composantes
               Prec=double(3),    # D x D x K : précisions des composantes
               pi=double(1),      # K : poids du mélange (somme ≈ 1)
               grid=double(2),    # M x 2 : grille 2D
               omega=double(1),   # M : poids de quadrature
               A=double(0),       # scalaire : facteur d’aire (dx*dy)
               K=integer(0),      # nombre de composantes
               d1=integer(0),     # indice de la 1re dimension
               d2=integer(0),     # indice de la 2e dimension
               M=integer(0)){     # nombre de points de grille
    returnType(double(0))
    
    # Buffers 2D
    mean2 <- numeric(2)            # moyenne projetée (d1, d2)
    Prec2 <- matrix(0.0, 2, 2)     # précision projetée 2x2
    x2    <- numeric(2)            # point de grille courant
    sumZ  <- 0.0                   # accumulation pour la quadrature
    
    for(m in 1:M){
      # point m de la grille
      x2[1] <- grid[m,1]; x2[2] <- grid[m,2]
      mix <- 0.0
      
      # mélange gaussien au point x2
      for(k in 1:K){
        # projection 2D de la moyenne et de la précision (sous-matrice)
        mean2[1] <- mu[k,d1]; mean2[2] <- mu[k,d2]
        Prec2[1,1] <- Prec[d1,d1,k]; Prec2[1,2] <- Prec[d1,d2,k]
        Prec2[2,1] <- Prec[d2,d1,k]; Prec2[2,2] <- Prec[d2,d2,k]
        
        # densité 2D de la composante k en x2 (écriture précision)
        mix <- mix + pi[k] * dmvnorm_nimble(x2, mean2, Prec2, FALSE)
      }
      
      # ajout quadrature pondérée
      sumZ <- sumZ + omega[m] * mix
    }
    
    # Échelle d’aire + protection contre l’underflow
    Z <- A * sumZ
    if (Z < 1e-300) Z <- 1e-300
    
    # renvoie log Z
    return(log(Z))
  }
)


# ======================= 6) Modèle corrigé par Ω (4D) =======================

code_corr4D <- nimbleCode({
  # ------------------ 1) Constantes de normalisation par saison (2D) ------------------
  # Été : intégration du mélange sur (d1, d2) = (1, 2) avec les poids d’été
  logZ_su <- logZ_season(mu[1:K,1:D], Prec[1:D,1:D,1:K], pi[1:K],
                         grid[1:M,1:2], omega_su[1:M], A, K, 1, 2, M)
  # Hiver : intégration du mélange sur (d1, d2) = (3, 4) avec les poids d’hiver
  logZ_wi <- logZ_season(mu[1:K,1:D], Prec[1:D,1:D,1:K], pi[1:K],
                         grid[1:M,1:2], omega_wi[1:M], A, K, 3, 4, M)
  logZ[1] <- logZ_su
  logZ[2] <- logZ_wi
  
  # ------------------ 2) Vraisemblance point par point (mélange gaussien 4D) ------------------
  
  for (i in 1:N) {
    for (k in 1:K) {
      dens[i,k] <- pi[k] * exp(dmvnorm_nimble(X[i,1:D], mu[k,1:D], Prec[1:D,1:D,k], TRUE))
    }
    
    # Densité totale du mélange en i (somme sur les composantes)
    mixdens[i] <- sum(dens[i,1:K])
    
    # Log-vraisemblance normalisée par la constante de la saison de l’observation i
    ll[i] <- log(Omega[i] * mixdens[i] + 1e-300) - logZ[sampling.season[i]]
    
    # ------------------ 3) Contrôle numérique et borne des probabilités ------------------
    # Stabilisation : recentrer par C_ub puis écrêter dans [p_min, p_max1]
    p_raw[i]     <- exp(ll[i] - C_ub)
    p_clip_hi[i] <- min(p_raw[i], p_max1)
    p[i]         <- max(p_clip_hi[i], p_min)
    ones[i] ~ dbern(p[i])
  }
  
  # ------------------ 4) Priors ------------------
  
  pi[1:K] ~ ddirch(alpha[1:K])
  
  for (k in 1:K) {
    mu[k,1:D] ~ dmnorm(mu0[1:D], prec = Prec0[1:D,1:D])
  }
  
  # Bloc-diagonal (précisions été/hiver séparées, sans covariation inter-saison)
  for (k in 1:K) {
    Prec_su[1:2,1:2,k] ~ dwish(R_su[1:2,1:2], df_su)
    Prec_wi[1:2,1:2,k] ~ dwish(R_wi[1:2,1:2], df_wi)
    
    for(i in 1:2) for(j in 1:2) Prec[i,   j,   k] <- Prec_su[i,j,k]
    for(i in 1:2) for(j in 1:2) Prec[i+2, j+2, k] <- Prec_wi[i,j,k]
    for(i in 1:2) for(j in 3:4) Prec[i,   j,   k] <- 0
    for(i in 3:4) for(j in 1:2) Prec[i,   j,   k] <- 0
  }
})

# ======================= 7) Données / constantes =======================

# Paramètres de stabilisation numérique et bornes des probabilités
C_ub   <- 200       # constante de recentrage pour éviter les dépassements numériques dans exp(ll - C_ub)
p_min  <- 1e-300    # probabilité minimale autorisée (évite p = 0 et log(0))
p_max1 <- 1 - 1e-12 # probabilité maximale autorisée (évite p = 1 exact, overflow ou log(1-p)=0)


# Constantes
constants_corr <- list(
  N=N, D=D, K=K,
  M=M, A=A_cell,
  alpha=rep(1,K),
  mu0=rep(0,D),
  Prec0=diag(D)*1e-2,
  R_su=diag(2), R_wi=diag(2),
  df_su=4, df_wi=4,
  # bornes du ones-trick
  p_min=p_min, p_max1=p_max1, C_ub=C_ub,
  sampling.season=as.integer(sampling.season)
)

# Données
data_corr <- list(
  X=X,
  Omega=Omega_vec,
  ones=rep(1L,N),
  grid=grid_mat,
  omega_su=omega_su_grid,
  omega_wi=omega_wi_grid
)

# ======================= 8) Initialisations =======================
# Barycentre global et écarts-types
center_global <- colMeans(X4_all)
sd_global <- apply(X4_all, 2, sd)

# Initialisation des moyennes : autour du barycentre (jitter)
mu_init <- matrix(NA_real_, nrow = K, ncol = D)
for (k in 1:K) {
  mu_init[k, ] <- center_global +
    c(rnorm(1, 0, 0.2 * sd_global[1]),
      rnorm(1, 0, 0.2 * sd_global[2]),
      rnorm(1, 0, 0.2 * sd_global[3]),
      rnorm(1, 0, 0.2 * sd_global[4]))
}

# Précisions initiales : inverse des covariances moyennes (bloquées été/hiver)
Sigma_su_mean <- Reduce("+", lapply(Sigma_true, function(S) S[1:2,1:2])) / K
Sigma_wi_mean <- Reduce("+", lapply(Sigma_true, function(S) S[3:4,3:4])) / K
Prec_su_init <- array(solve(Sigma_su_mean), dim = c(2, 2, K))
Prec_wi_init <- array(solve(Sigma_wi_mean), dim = c(2, 2, K))

# Léger jitter sur les précisions
for (k in 1:K) {
  Prec_su_init[,,k] <- Prec_su_init[,,k] + diag(runif(2, 0, 0.05))
  Prec_wi_init[,,k] <- Prec_wi_init[,,k] + diag(runif(2, 0, 0.05))
}

# Liste finale d’initialisation
inits1 <- list(
  mu = mu_init,
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
samples_naif4 <- runMCMC(cmcmc_n4, niter=3000, nburnin=1000, thin=5, setSeed=TRUE)

cat("[Naïf 4D] Terminé.\n")
