


set.seed(123)

suppressPackageStartupMessages({
  library(MASS)
  library(Matrix)   # nearPD
  library(terra)
  library(ggplot2)
  library(patchwork)
  library(nimble)
  library(mvtnorm)
  library(viridisLite)
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
N_tot <- 5000
pi_true <- c(0.20, 0.45, 0.35)

# Moyennes 4D (lon_été, lat_été, lon_hiver, lat_hiver) — stéréotypées
mu <- list(
  c(-3.5,  0.0,  -3.0,  1.2),  # C1 : été = gauche (vertical), hiver = gauche (diag +45°)
  c( 0.0, -0.6,   0.2, -0.4),  # C2 : été = diagonal -45° (centre), hiver = ~vertical centre
  c( 3.5,  0.0,   3.2,  0.3)   # C3 : été = droite (vertical), hiver = droite (horizontal)
)
# Décalage global hiver (axe y) pour séparer visuellement les saisons
for(k in 1:K) mu[[k]][4] <- mu[[k]][4] - 30

# Covariances 4x4 complètes (corrélations été↔hiver) — stéréotypées
## C1 : Été vertical (sd_x petit, sd_y grand), Hiver diagonal +45°
sd1 <- c(0.7, 2.8, 2.6, 1.9)  # (x_su, y_su, x_wi, y_wi)
R1 <- matrix(c(
  1.00,  0.00,  0.55, -0.15,
  0.00,  1.00,  0.10,  0.30,
  0.55,  0.10,  1.00,  0.65,
  -0.15,  0.30,  0.65,  1.00
), 4, 4, byrow = TRUE)
Sigma1 <- build_sigma(sd1, R1)

## C2 : Été diagonal -45° (corr < 0), Hiver vertical (sd_x petit, sd_y grand)
sd2 <- c(2.6, 2.4, 0.8, 3.0)
R2 <- matrix(c(
  1.00, -0.70,  0.40, -0.35,
  -0.70,  1.00, -0.10,  0.45,
  0.40, -0.10,  1.00,  0.00,
  -0.35,  0.45,  0.00,  1.00
), 4, 4, byrow = TRUE)
Sigma2 <- build_sigma(sd2, R2)

## C3 : Été vertical (à droite), Hiver horizontal (sd_x grand, sd_y petit)
sd3 <- c(0.8, 3.2, 3.1, 0.9)
R3 <- matrix(c(
  1.00,  0.00,  0.50,  0.05,
  0.00,  1.00, -0.20,  0.25,
  0.50, -0.20,  1.00,  0.05,
  0.05,  0.25,  0.05,  1.00
), 4, 4, byrow = TRUE)
Sigma3 <- build_sigma(sd3, R3)

Sigma <- list(Sigma1, Sigma2, Sigma3)
true_Sigma <- Sigma  # pour les plots

# Simulation 4D (sans Ω)
z <- sample(1:K, N_tot, replace = TRUE, prob = pi_true)
X4_all <- t(vapply(
  1:N_tot, function(i) MASS::mvrnorm(1, mu = mu[[z[i]]], Sigma = Sigma[[z[i]]]),
  numeric(4)
))
colnames(X4_all) <- c("x_su","y_su","x_wi","y_wi")

# Séparations 2D pour Ω/plots
X_su_all <- X4_all[, 1:2, drop = FALSE]
X_wi_all <- X4_all[, 3:4, drop = FALSE]

# ==================================================
# 2) Ω(x, saison) + Filtrage
# ==================================================
# ÉTÉ : identique à l’original (3 bosses)
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

# HIVER : 1 large tache explicite proche du cluster hiver C1 (décalée pour limiter le recouvrement)
build_omega_winter <- function(xmin, xmax, ymin, ymax, nx = 80, ny = 80,
                               mu, enlarge = 1.0){
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

# Emprise commune
pad <- 1
xmin_all <- min(c(X_su_all[,1], X_wi_all[,1])) - pad
xmax_all <- max(c(X_su_all[,1], X_wi_all[,1])) + pad
ymin_all <- min(c(X_su_all[,2], X_wi_all[,2])) - pad
ymax_all <- max(c(X_su_all[,2], X_wi_all[,2])) + pad

# Rasters Ω
r_su <- build_omega_summer(xmin_all, xmax_all, ymin_all, ymax_all, nx = 80, ny = 80)
r_wi <- build_omega_winter(xmin_all, xmax_all, ymin_all, ymax_all, nx = 80, ny = 80, mu = mu)


Omega_su_all <- terra::extract(r_su, X_su_all)[,1]
Omega_wi_all <- terra::extract(r_wi, X_wi_all)[,1]

# ==================================================
# 3) Échantillonnage Été / Hiver
# ==================================================
sel_ete <- runif(N_tot) < Omega_su_all
sel_hiver <- runif(N_tot) < Omega_wi_all
obs_ete <- data.frame(id=which(sel_ete), saison="été",
                      x_su=X_su_all[sel_ete,1], y_su=X_su_all[sel_ete,2],
                      x_wi=X_wi_all[sel_ete,1], y_wi=X_wi_all[sel_ete,2],
                      Omega=Omega_su_all[sel_ete])
obs_hiver <- data.frame(id=which(sel_hiver), saison="hiver",
                        x_su=X_su_all[sel_hiver,1], y_su=X_su_all[sel_hiver,2],
                        x_wi=X_wi_all[sel_hiver,1], y_wi=X_wi_all[sel_hiver,2],
                        Omega=Omega_wi_all[sel_hiver])
obs_long <- rbind(obs_ete, obs_hiver)
N <- nrow(obs_long)

# ==================================================
# 4) NIMBLE modèle corrigé (zeros-trick stabilisé)
# ==================================================

# Calcul la densité normale multivariée à une point précis pour une GMM quelconque
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
    phi[i] <- pmax(-ll_clip[i] + 1e-6, 1e-6)
  }
  
  pi[1:K] ~ ddirch(alpha[1:K])
  for(k in 1:K){
    mu[k,1:D] ~ dmnorm(mu0[1:D], prec=Prec0[1:D,1:D])
    Prec[1:D,1:D,k] ~ dwish(R[1:D,1:D], df)
  }
})




# ==================================================
# 5) Données pour NIMBLE
# ==================================================
D <- 4L
X <- as.matrix(obs_long[,c("x_su","y_su","x_wi","y_wi")])
Omega <- obs_long$Omega
ss <- ifelse(obs_long$saison=="été",1L,2L)
zeros <- rep(0L,N)

su_df <- as.data.frame(r_su, xy=TRUE); names(su_df)[3] <- "omega"
wi_df <- as.data.frame(r_wi, xy=TRUE); names(wi_df)[3] <- "omega"
A_su <- prod(res(r_su)); A_wi <- prod(res(r_wi))
M_su <- 4000L; M_wi <- 4000L
grid_su <- as.matrix(su_df[sample.int(nrow(su_df), M_su, TRUE), c("x","y")])
grid_wi <- as.matrix(wi_df[sample.int(nrow(wi_df), M_wi, TRUE), c("x","y")])
omega_su_grid <- su_df$omega[sample.int(nrow(su_df), M_su, TRUE)]
omega_wi_grid <- wi_df$omega[sample.int(nrow(wi_df), M_wi, TRUE)]

constants <- list(N=N, D=D, K=K, M_su=M_su, M_wi=M_wi,
                  A_su=A_su, A_wi=A_wi,
                  alpha=rep(1,K), mu0=rep(0,D),
                  Prec0=diag(D)*1e-2, R=diag(D), df=D+2)
data_list <- list(X=X, Omega=Omega, ss=ss, zeros=zeros,
                  grid_su=grid_su, grid_wi=grid_wi,
                  omega_su_grid=omega_su_grid, omega_wi_grid=omega_wi_grid)

# ==================================================
# 6) Inits barycentre + bruit
# ==================================================
barX <- colMeans(X)
mu_init <- matrix(NA,K,D)
for(k in 1:K) mu_init[k,] <- barX + rnorm(D,0,1)
Prec_init <- array(0,dim=c(D,D,K))
for(k in 1:K) Prec_init[,,k] <- diag(D)*0.5
inits1 <- list(mu=mu_init,Prec=Prec_init,pi=rep(1/K,K))
inits2 <- list(mu=mu_init+rnorm(K*D,0,0.5),Prec=Prec_init,pi=rep(1/K,K))

# ==================================================
# 7) Compilation et exécution (2 chaînes)
# ==================================================
# 7) Compilation et exécution du MCMC (2 chaînes)
# ==================================================
model <- nimbleModel(code_corrige,
                     data = data_list,
                     constants = constants,
                     inits = inits1,
                     check = FALSE)

Sys.setenv("CXXFLAGS"="-Wno-format -Wno-error")
cmodel <- compileNimble(model, showCompilerOutput=TRUE)
conf <- configureMCMC(model, monitors = c("mu", "pi", "Prec"))
mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc, project = model)

cat("➡️ Lancement de la chaîne 1...\n")
samples1 <- runMCMC(cmcmc, niter = 2000, nburnin = 500, thin = 2,
                    inits = inits1, setSeed = 100)

cat("➡️ Lancement de la chaîne 2...\n")
samples2 <- runMCMC(cmcmc, niter = 2000, nburnin = 500, thin = 2,
                    inits = inits2, setSeed = 200)



























