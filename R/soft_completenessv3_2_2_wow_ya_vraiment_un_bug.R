# ==================================================
# Modèle NIMBLE (GMM 4D + zeros-trick + Ω)
#  Solution A : pas de corrélations été↔hiver (bloc-diagonal)
#  + prior ordonné sur mu[,1] (labels stables)
#  + zeros-trick robuste (phi >= 0)
# ==================================================


set.seed(123)
suppressPackageStartupMessages({
  library(MASS)
  library(Matrix)
  library(terra)
  library(nimble)
  library(mclust)
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
pi_true <- c(0.20, 0.35, 0.45)

mu <- list(
  c(-3.5,  0.0,  0.5 - 3.0,  0.0),
  c( 0.0, -0.6,  1.0 - 3.0,  10),
  c( 3.5,  0.0,  6.0 - 3.0,  3.0)
)
for(k in 1:K) mu[[k]][4] <- mu[[k]][4] - 30

# Bloc-diagonal : (1,2) ⟂ (3,4)
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
build_omega_summer <- function(xmin, xmax, ymin, ymax, nx = 60, ny = 60, global_intensity = 0.5){
  r <- rast(nrows = ny, ncols = nx, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
  xy <- as.data.frame(xyFromCell(r, 1:ncell(r))); colnames(xy) <- c("x","y")
  omega_fun <- function(x, y){
    base <- 0.04
    bumpL <- 0.90 * gauss2d(x, y, c(-3.8, 2.6), diag(c(0.60, 0.60)))
    bumpR <- 0.90 * gauss2d(x, y, c( 3.8, 2.6), diag(c(0.65, 0.65)))
    bumpC <- 1.00 * gauss2d(x, y, c( 0.0,-3.2), diag(c(6.0, 0.9)))
    pmax(0, pmin(base + bumpL + bumpR + bumpC, 0.70))
  }
  r[] <- with(xy, omega_fun(x, y)) * global_intensity
  r
}
build_omega_winter <- function(xmin, xmax, ymin, ymax, nx = 60, ny = 60, mu, enlarge = 2.0, global_intensity = 0.3){
  center_wi_1 <- c(0, -26)
  Sigma_wi_1  <- make_cov2D(angle_deg =  1, var_major = 2.0*enlarge, var_minor = 2.0*enlarge)
  base_wi <- 0.04; cap_wi <- 0.90; amp1 <- 1.00
  r <- rast(nrows = ny, ncols = nx, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
  xy <- as.data.frame(xyFromCell(r, 1:ncell(r))); colnames(xy) <- c("x","y")
  omega_fun <- function(x, y){
    val <- base_wi + amp1 * gauss2d(x, y, center_wi_1, Sigma_wi_1)
    pmax(0, pmin(val, cap_wi))
  }
  r[] <- with(xy, omega_fun(x, y)) * global_intensity
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
Omega_su_all[!is.finite(Omega_su_all)] <- 0
Omega_wi_all[!is.finite(Omega_wi_all)] <- 0

sel_ete   <- runif(N_tot) < Omega_su_all
sel_hiver <- runif(N_tot) < Omega_wi_all
obs_ete   <- data.frame(id=which(sel_ete), saison="été",
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
# 3) Grille LxL appariée (intégration)
# ==================================================
D <- 4L
X <- as.matrix(obs_long[,c("x_su","y_su","x_wi","y_wi")])
Omega <- obs_long$Omega
ss <- ifelse(obs_long$saison=="été",1L,2L)
zeros <- rep(0L,N)

L <- 60L
ext_all <- ext(r_su)
xmin <- ext_all[1]; xmax <- ext_all[2]
ymin <- ext_all[3]; ymax <- ext_all[4]
dx <- (xmax - xmin) / L
dy <- (ymax - ymin) / L
A_cell <- dx * dy
xs <- xmin + (0:(L-1) + 0.5) * dx
ys <- ymin + (0:(L-1) + 0.5) * dy
grid_mat <- as.matrix(expand.grid(x = xs, y = ys))

omega_su_grid <- as.numeric(terra::extract(r_su, grid_mat, method = "bilinear")[,1])
omega_wi_grid <- as.numeric(terra::extract(r_wi, grid_mat, method = "bilinear")[,1])
omega_su_grid[!is.finite(omega_su_grid)] <- 0
omega_wi_grid[!is.finite(omega_wi_grid)] <- 0

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
    Z <- sumZ * A
    if(!(Z > 0.0)) Z <- 1e-300      # blindage anti-NA/<=0
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
      dens[i,k] <- pi[k] * exp(dmvnorm_nimble(X[i,1:D], mu[k,1:D], Prec[1:D,1:D,k], TRUE))
    }
    mixdens[i] <- sum(dens[i,1:K])
    ll[i]      <- log(Omega[i]) + log(mixdens[i]) - logZ[ss[i]] 
    
    # zeros-trick robuste : phi >= 0
    phi[i]     <- (Cclip - min(ll[i], Cclip)) + 1e-8
    zeros[i]   ~ dpois(phi[i])
  }
  
  ll_global <- sum(ll[1:N])
  
  ll_min <-min(ll[1:N])
  ll_mean <-mean(ll[1:N])
  ll_max <-max(ll[1:N])
  
  log_mixdens_min <- min(log(mixdens[1:N]))
  log_mixdens_mean <- mean(log(mixdens[1:N]))
  log_mixdens_max <- max(log(mixdens[1:N]))
  
  log_Omega_min <- min(log(Omega[1:N]))
  log_Omega_mean <- mean(log(Omega[1:N]))
  log_Omegas_max <- max(log(Omega[1:N]))
  
  phi_min <- min(phi[1:N])
  phi_mean <- mean(phi[1:N])
  phi_max <- max(phi[1:N])
  
  pi_raw[1:K] ~ ddirch(alpha[1:K])
  for (k in 1:K) {
    pi[k] <- (1 - K * 0.15) * pi_raw[k] + 0.15
  }
  
  ## Prior ordonné sur mu[,1]
  mu1_base ~ dnorm(mu0[1], tau = prec_mu1_base)
  mu[1,1]  <- mu1_base
  for (k in 1:(K-1)) {
    gap[k]  ~ dexp(lambda_gap)
    mu[k+1,1] <- mu[k,1] + gap[k]
  }
  
  ## Priors classiques autres dimensions
  for (k in 1:K) for (d in 2:D) {
    mu[k,d] ~ dnorm(mu0[d], tau = Prec0[d,d])
  }
  
  for (k in 1:K) {
    Prec_su[1:2,1:2,k] ~ dwish(R_su[1:2,1:2], df_su)
    Prec_wi[1:2,1:2,k] ~ dwish(R_wi[1:2,1:2], df_wi)
    
    for(i in 1:2) for(j in 1:2) Prec[i,   j,   k] <- Prec_su[i,j,k]
    for(i in 1:2) for(j in 1:2) Prec[i+2, j+2, k] <- Prec_wi[i,j,k]
    for(i in 1:2) for(j in 3:4) Prec[i,   j,   k] <- 0
    for(i in 3:4) for(j in 1:2) Prec[i,   j,   k] <- 0
  }
})

code_naif <- nimbleCode({
  for(i in 1:N){
    for(k in 1:K){
      dens[i,k] <- pi[k] * exp(dmvnorm_nimble(X[i,1:D], mu[k,1:D],
                                              Prec[1:D,1:D,k], TRUE))
    }
    mixdens[i] <- sum(dens[i,1:K])
    ll[i]    <- log(mixdens[i] + 1e-300)
    
    # même zeros-trick robuste pour éviter phi<0
    phi[i]   <- (Cclip - min(ll[i], Cclip)) + 1e-8
    zeros[i] ~ dpois(phi[i])
  }
  
  ll_global <- sum(ll[1:N])
  
  
  pi[1:K] ~ ddirch(alpha[1:K])
  
  for (k in 1:K) for (d in 1:D) {
    mu[k,d] ~ dnorm(mu0[d], tau = Prec0[d,d])
  }
  
  for (k in 1:K) {
    Prec_su[1:2,1:2,k] ~ dwish(R_su[1:2,1:2], df_su)
    Prec_wi[1:2,1:2,k] ~ dwish(R_wi[1:2,1:2], df_wi)
    
    for(i in 1:2) for(j in 1:2) Prec[i,   j,   k] <- Prec_su[i,j,k]
    for(i in 1:2) for(j in 1:2) Prec[i+2, j+2, k] <- Prec_wi[i,j,k]
    for(i in 1:2) for(j in 3:4) Prec[i,   j,   k] <- 0
    for(i in 3:4) for(j in 1:2) Prec[i,   j,   k] <- 0
  }
})

# ==================================================
# 5) Données/constantes pour NIMBLE
# ==================================================
constants <- list(
  N=N, D=D, K=K,
  M_su=M_su, M_wi=M_wi,
  A_su=A_su, A_wi=A_wi,
  alpha=rep(1,K), mu0=rep(0,D),
  Prec0=diag(D)*1e-2,
  R_su=diag(2), R_wi=diag(2),
  df_su=4, df_wi=4,
  prec_mu1_base = 1/25,
  lambda_gap    = 1.0,
  Cclip = 700,                 # utilisé dans les deux modèles
  ss=as.integer(ss)
)

data_list <- list(
  X=X,
  Omega=Omega,
  zeros=as.integer(zeros),
  grid_su=grid_su, grid_wi=grid_wi,
  omega_su_grid=omega_su_grid, omega_wi_grid=omega_wi_grid
)

constants_naif <- list(
  N = N, D = D, K = K,
  alpha = rep(1, K),
  mu0   = rep(0, D),
  Prec0 = diag(D) * 1e-2,
  R_su = diag(2), R_wi = diag(2),
  df_su = 4, df_wi = 4,
  Cclip = 700                  # mêmes conventions
)

data_list_naif <- list(
  X = X,
  zeros = as.integer(zeros)
)

# ==================================================
# 6) Inits
# ==================================================
em <- Mclust(X, G = K, modelNames = "VVV", verbose = FALSE)
mu_init <- t(em$parameters$mean)

pi_init <- as.numeric(em$parameters$pro)
pi_init <- pmax(pi_init, 1e-6); pi_init <- pi_init / sum(pi_init)

Prec_init <- array(NA_real_, dim = c(D, D, K))
for (k in 1:K) {
  Sk <- em$parameters$variance$sigma[,,k]
  Sk <- (Sk + t(Sk)) / 2
  eig <- eigen(Sk, symmetric = TRUE, only.values = TRUE)$values
  if (any(!is.finite(eig)) || any(eig <= 1e-8)) Sk <- Sk + diag(1e-4, D)
  Qk <- solve(Sk)
  Qk[1:2,3:4] <- 0; Qk[3:4,1:2] <- 0   # bloc-diagonal forcé pour l'init
  Prec_init[,,k] <- Qk
}

# Ordonner les composantes par mu[,1]
ord <- order(mu_init[,1])
mu_init <- mu_init[ord, , drop = FALSE]
Prec_init <- Prec_init[,, ord, drop = FALSE]
pi_init <- pi_init[ord]

# Inits pour les nœuds stochastiques Prec_su / Prec_wi
Prec_su_init <- array(NA_real_, dim = c(2,2,K))
Prec_wi_init <- array(NA_real_, dim = c(2,2,K))
for(k in 1:K){
  Prec_su_init[,,k] <- Prec_init[1:2,1:2,k]
  Prec_wi_init[,,k] <- Prec_init[3:4,3:4,k]
}

# Inits pour le prior ordonné (mu1_base, gaps>0)
mu_sd <- 0.15
mu1_base_init <- mu_init[1,1]
gap_init <- pmax(1e-3, diff(mu_init[,1]))

inits1 <- list(
  mu   = mu_init + matrix(rnorm(K*D, 0, mu_sd), K, D),
  pi   = pi_init,
  Prec_su = Prec_su_init,
  Prec_wi = Prec_wi_init,
  mu1_base = mu1_base_init,
  gap = gap_init
)
pi_sd <- 0.02
pi2 <- pmax(1e-6, pi_init + rnorm(K, 0, pi_sd)); pi2 <- pi2/sum(pi2)
inits2 <- list(
  mu   = mu_init + matrix(rnorm(K*D, 0, mu_sd), K, D),
  pi   = pi2,
  Prec_su = Prec_su_init,
  Prec_wi = Prec_wi_init,
  mu1_base = mu1_base_init + rnorm(1, 0, 0.05),
  gap = pmax(1e-3, gap_init + rnorm(K-1, 0, 0.05))
)

# ==================================================
# 7) Modèle + compile + MCMC (2 chaînes)
# ==================================================
model <- nimbleModel(code_corrige, data = data_list, constants = constants,
                     inits = inits1, check = FALSE)

conf  <- configureMCMC(model, monitors = c("mu","pi","Prec", "ll_global", "logZ",
                                           "ll_min",
                                           "ll_mean",
                                           "ll_max",
                                           
                                           "log_mixdens_min",
                                           "log_mixdens_mean",
                                           "log_mixdens_max",
                                           
                                           "log_Omega_min",
                                           "log_Omega_mean",
                                           "log_Omegas_max",
                                           
                                           "phi_min",
                                           "phi_mean",
                                           "phi_max"
                                           
))
mcmc  <- buildMCMC(conf)

Cmodel <- nimble::compileNimble(model, showCompilerOutput = FALSE)
cmcmc  <- nimble::compileNimble(mcmc, project = Cmodel, showCompilerOutput = FALSE)

cat("➡️ Chaîne 1...\n")
samples1 <- runMCMC(cmcmc, niter = 1200, nburnin = 2, thin = 2,
                    inits = inits1, setSeed = 100)

cat("➡️ Chaîne 2...\n")
samples2 <- runMCMC(cmcmc, niter = 1200, nburnin = 2, thin = 2,
                    inits = inits2, setSeed = 200)

# =============================
# Modèle NAÏF (même format blocs)
# =============================
model_naif <- nimbleModel(code_naif,
                          data      = data_list_naif,
                          constants = constants_naif,
                          inits     = inits1, check = FALSE)

conf_naif  <- configureMCMC(model_naif, monitors = c("mu","pi","Prec", "ll_global"))
mcmc_naif  <- buildMCMC(conf_naif)

Cmodel_naif <- nimble::compileNimble(model_naif, showCompilerOutput = FALSE)
cmcmc_naif  <- nimble::compileNimble(mcmc_naif, project = Cmodel_naif, showCompilerOutput = FALSE)

cat("➡️ Chaîne 1 (naïf)...\n")
samples_naif1 <- runMCMC(cmcmc_naif, niter = 1000, nburnin = 2, thin = 2,
                         inits = inits1, setSeed = 101)

cat("➡️ Chaîne 2 (naïf)...\n")
samples_naif2 <- runMCMC(cmcmc_naif, niter = 1000, nburnin = 2, thin = 2,
                         inits = inits2, setSeed = 202)
