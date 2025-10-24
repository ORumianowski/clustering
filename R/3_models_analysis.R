# =============================================================
# Script 3/4 — Analyse par modèles (NIMBLE)
# =============================================================
# Réutilise:
#   - Script 1: contextes + Ω + helpers
#   - Script 2: simulation (datasets)
# Contenu:
#   1) Paramètres de tuning (quadrature, MCMC)
#   2) Fonctions NIMBLE utilitaires
#   3) Modèles: naïf 4D ; corrigé Ω 4D
#   4) Préparation des entrées pour chaque modèle
#   5) Exécuteurs (2 chaînes par défaut)
#   6) Exemples d'utilisation (commentés)
# =============================================================

suppressPackageStartupMessages({
  library(nimble)
  library(MASS)
  library(terra)
})

# --- IMPORTANT: réutiliser Scripts 1 & 2 si pas déjà chargés ---
if (!exists("get_context_params")) {
  source("R/1_contexts_omegas.R")
}
if (!exists("simulate_dataset")) {
  source("R/2_simulate_datasets.R")
}

# =============================================================
# 1) Paramètres de tuning (par défaut)
# =============================================================

gmm_tuning <- list(
  # Quadrature pour intégrales de normalisation saisonnières (2D)
  L_grid = 60L,        # nombre de points par côté
  # Bornes stabilisation ones-trick
  C_ub   = 200,        # recentrage pour exp(ll - C_ub)
  p_min  = 1e-300,
  p_max1 = 1 - 1e-12,
  # MCMC
  niter   = 3000,
  nburnin = 1000,
  thin    = 5
)

# =============================================================
# 2) Fonctions NIMBLE utilitaires
# =============================================================

dmvnorm_nimble <- nimbleFunction(
  run = function(x=double(1), mean=double(1), Prec=double(2), log=logical(0, default=FALSE)){
    returnType(double(0))
    Dloc <- length(x)
    xm <- x - mean
    qf <- inprod(xm, Prec %*% xm)
    U <- chol(Prec)
    ldet <- 2 * sum(log(diag(U)))
    logdens <- 0.5 * ldet - 0.5 * Dloc * log(2*pi) - 0.5 * qf
    if (log) return(logdens) else return(exp(logdens))
  }
)

logZ_season <- nimbleFunction(
  run=function(mu=double(2), Prec=double(3), pi=double(1), grid=double(2), omega=double(1), A=double(0),
               K=integer(0), d1=integer(0), d2=integer(0), M=integer(0)){
    returnType(double(0))
    mean2 <- numeric(2)
    Prec2 <- matrix(0.0, 2, 2)
    x2    <- numeric(2)
    sumZ  <- 0.0
    for(m in 1:M){
      x2[1] <- grid[m,1]; x2[2] <- grid[m,2]
      mix <- 0.0
      for(k in 1:K){
        mean2[1] <- mu[k,d1]; mean2[2] <- mu[k,d2]
        Prec2[1,1] <- Prec[d1,d1,k]; Prec2[1,2] <- Prec[d1,d2,k]
        Prec2[2,1] <- Prec[d2,d1,k]; Prec2[2,2] <- Prec[d2,d2,k]
        mix <- mix + pi[k] * dmvnorm_nimble(x2, mean2, Prec2, FALSE)
      }
      sumZ <- sumZ + omega[m] * mix
    }
    Z <- A * sumZ
    if (Z < 1e-300) Z <- 1e-300
    return(log(Z))
  }
)

# =============================================================
# 3) Modèles NIMBLE
# =============================================================

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

code_corr4D <- nimbleCode({
  # Constantes de normalisation par saison
  logZ_su <- logZ_season(mu[1:K,1:D], Prec[1:D,1:D,1:K], pi[1:K], grid[1:M,1:2], omega_su[1:M], A, K, 1, 2, M)
  logZ_wi <- logZ_season(mu[1:K,1:D], Prec[1:D,1:D,1:K], pi[1:K], grid[1:M,1:2], omega_wi[1:M], A, K, 3, 4, M)
  logZ[1] <- logZ_su
  logZ[2] <- logZ_wi
  
  for (i in 1:N) {
    for (k in 1:K) {
      dens[i,k] <- pi[k] * exp(dmvnorm_nimble(X[i,1:D], mu[k,1:D], Prec[1:D,1:D,k], TRUE))
    }
    mixdens[i] <- sum(dens[i,1:K])
    ll[i] <- log(Omega[i] * mixdens[i] + 1e-300) - logZ[sampling.season[i]]
    p_raw[i]     <- exp(ll[i] - C_ub)
    p_clip_hi[i] <- min(p_raw[i], p_max1)
    p[i]         <- max(p_clip_hi[i], p_min)
    ones[i] ~ dbern(p[i])
  }
  
  pi[1:K] ~ ddirch(alpha[1:K])
  for (k in 1:K) {
    mu[k,1:D] ~ dmnorm(mu0[1:D], prec = Prec0[1:D,1:D])
    # Bloc-diagonal été/hiver
    Prec_su[1:2,1:2,k] ~ dwish(R_su[1:2,1:2], df_su)
    Prec_wi[1:2,1:2,k] ~ dwish(R_wi[1:2,1:2], df_wi)
    for(i in 1:2) for(j in 1:2) Prec[i,   j,   k] <- Prec_su[i,j,k]
    for(i in 1:2) for(j in 1:2) Prec[i+2, j+2, k] <- Prec_wi[i,j,k]
    for(i in 1:2) for(j in 3:4) Prec[i,   j,   k] <- 0
    for(i in 3:4) for(j in 1:2) Prec[i,   j,   k] <- 0
  }
})

# =============================================================
# 4) Préparation des entrées
# =============================================================

make_grid_quadrature <- function(r_ref, L=60){
  ext_all <- ext(r_ref)
  xmin <- ext_all[1]; xmax <- ext_all[2]
  ymin <- ext_all[3]; ymax <- ext_all[4]
  dx <- (xmax - xmin) / L; dy <- (ymax - ymin) / L
  xs <- xmin + (0:(L-1) + 0.5) * dx
  ys <- ymin + (0:(L-1) + 0.5) * dy
  grid_mat <- as.matrix(expand.grid(x = xs, y = ys))
  list(grid=grid_mat, A_cell=dx*dy)
}

# ---- Entrées pour le modèle corrigé (Ω utilisé pour l'analyse configurable) ----
prepare_inputs_corr <- function(sim, analysis_omega_version=c("original","deg1","deg2"), L=gmm_tuning$L_grid){
  analysis_omega_version <- match.arg(analysis_omega_version)
  
  # Ω pour l'analyse (peut être dégradé)
  om_analysis <- make_omega_pair(analysis_omega_version, n_pix_side = ncol(sim$r_su))
  r_su_use <- om_analysis$summer
  r_wi_use <- om_analysis$winter
  
  df <- sim$data_long
  X <- as.matrix(df[,c("x_su","y_su","x_wi","y_wi")])
  sampling.season <- ifelse(df$saison=="ete", 1L, 2L)
  Omega_vec <- df$Omega
  N <- nrow(df); D <- 4L; K <- 3L
  
  g <- make_grid_quadrature(r_su_use, L=L)
  grid_mat <- g$grid; A_cell <- g$A_cell; M <- nrow(grid_mat)
  
  omega_su_grid <- as.numeric(terra::extract(r_su_use, grid_mat, method = "bilinear")[,1])
  omega_wi_grid <- as.numeric(terra::extract(r_wi_use, grid_mat, method = "bilinear")[,1])
  
  constants <- list(
    N=N, D=D, K=K,
    M=M, A=A_cell,
    alpha=rep(1,K),
    mu0=rep(0,D),
    Prec0=diag(D)*1e-2,
    R_su=diag(2), R_wi=diag(2),
    df_su=4, df_wi=4,
    p_min=gmm_tuning$p_min, p_max1=gmm_tuning$p_max1, C_ub=gmm_tuning$C_ub,
    sampling.season=as.integer(sampling.season)
  )
  data <- list(
    X=X,
    Omega=Omega_vec,
    ones=rep(1L,N),
    grid=grid_mat,
    omega_su=omega_su_grid,
    omega_wi=omega_wi_grid
  )
  
  # Inits (autour de la data simulée)
  center_global <- colMeans(sim$X4_all)
  sd_global <- apply(sim$X4_all, 2, sd)
  mu_init <- matrix(NA_real_, nrow = K, ncol = D)
  for (k in 1:K) {
    mu_init[k, ] <- center_global + c(rnorm(1,0,0.2*sd_global[1]), rnorm(1,0,0.2*sd_global[2]),
                                      rnorm(1,0,0.2*sd_global[3]), rnorm(1,0,0.2*sd_global[4]))
  }
  Sigma_su_mean <- Reduce("+", lapply(sim$Sigma_true, function(S) S[1:2,1:2])) / length(sim$Sigma_true)
  Sigma_wi_mean <- Reduce("+", lapply(sim$Sigma_true, function(S) S[3:4,3:4])) / length(sim$Sigma_true)
  Prec_su_init <- array(solve(Sigma_su_mean), dim = c(2,2,K))
  Prec_wi_init <- array(solve(Sigma_wi_mean), dim = c(2,2,K))
  for (k in 1:K) {
    Prec_su_init[,,k] <- Prec_su_init[,,k] + diag(runif(2, 0, 0.05))
    Prec_wi_init[,,k] <- Prec_wi_init[,,k] + diag(runif(2, 0, 0.05))
  }
  inits <- list(mu = mu_init, pi = rep(1/K, K), Prec_su = Prec_su_init, Prec_wi = Prec_wi_init)
  
  list(constants=constants, data=data, inits=inits)
}

# ---- Entrées pour le modèle naïf ----
prepare_inputs_naif <- function(sim){
  df <- sim$data_long
  X <- as.matrix(df[,c("x_su","y_su","x_wi","y_wi")])
  N <- nrow(df); D <- 4L; K <- 3L
  
  center_global <- colMeans(sim$X4_all)
  sd_global <- apply(sim$X4_all, 2, sd)
  mu_init <- matrix(NA_real_, nrow = K, ncol = D)
  for (k in 1:K) {
    mu_init[k, ] <- center_global + c(rnorm(1,0,0.2*sd_global[1]), rnorm(1,0,0.2*sd_global[2]),
                                      rnorm(1,0,0.2*sd_global[3]), rnorm(1,0,0.2*sd_global[4]))
  }
  constants <- list(N=N, D=D, K=K, alpha=rep(1,K), mu0=rep(0,D), Prec0=diag(D)*1e-2, R4=diag(D), df4=D+2)
  data <- list(X=X)
  inits <- list(z=sample(1:K,N,TRUE), mu = mu_init, Prec = array(0, dim=c(D,D,K)), pi=rep(1/K,K))
  for(k in 1:K){ inits$Prec[,,k] <- diag(D) }
  list(constants=constants, data=data, inits=inits)
}

# =============================================================
# 5) Exécuteurs (2 chaînes)
# =============================================================

run_model_naif <- function(sim, niter=gmm_tuning$niter, nburnin=gmm_tuning$nburnin, thin=gmm_tuning$thin, seed_base=11){
  inp <- prepare_inputs_naif(sim)
  out <- vector("list", 2)
  for (ch in 1:2){
    set.seed(seed_base + ch)
    model <- nimbleModel(code_naif4D, data=inp$data, constants=inp$constants, inits=inp$inits, check=FALSE)
    cmodel <- compileNimble(model)
    conf <- configureMCMC(model, monitors=c("mu","pi","Prec"))
    mcmc <- buildMCMC(conf)
    cmcmc <- compileNimble(mcmc, project=cmodel)
    samples <- runMCMC(cmcmc, niter=niter, nburnin=nburnin, thin=thin, setSeed=TRUE)
    out[[ch]] <- samples
  }
  out
}

run_model_corr <- function(sim, analysis_omega_version=c("original","deg1","deg2"), L=gmm_tuning$L_grid,
                           niter=gmm_tuning$niter, nburnin=gmm_tuning$nburnin, thin=gmm_tuning$thin, seed_base=21){
  inp <- prepare_inputs_corr(sim, analysis_omega_version=analysis_omega_version, L=L)
  out <- vector("list", 2)
  for (ch in 1:2){
    set.seed(seed_base + ch)
    model <- nimbleModel(code_corr4D, data=inp$data, constants=inp$constants, inits=inp$inits, check=FALSE)
    cmodel <- compileNimble(model)
    conf <- configureMCMC(model, monitors=c("mu","pi","Prec","logZ"))
    mcmc <- buildMCMC(conf)
    cmcmc <- compileNimble(mcmc, project=cmodel)
    samples <- runMCMC(cmcmc, niter=niter, nburnin=nburnin, thin=thin, setSeed=TRUE)
    out[[ch]] <- samples
  }
  out
}

# =============================================================
# 6) Exemples d'utilisation (commentés)
# =============================================================

# ## 0) Charger scripts 1 & 2 si besoin
# source("1_contexts_omegas.R")
# source("2_simulate_datasets.R")
# 
# ## 1) Simuler (toujours avec Ω original) — choisir le contexte
# simC2 <- simulate_dataset("C2", N_tot=4000, seed=100)
# 
# ## 2) Lancer modèle naïf (2 chaînes)
# out_naif <- run_model_naif(simC2)
# 
# ## 3) Lancer modèle corrigé avec Ω utilisé pour l'analyse :
# ##    - "original" (oracle), "deg1" ou "deg2" (dégradés)
# out_corr_deg1 <- run_model_corr(simC2, analysis_omega_version="deg1", L=60)
# out_corr_deg2 <- run_model_corr(simC2, analysis_omega_version="deg2", L=60)
# out_corr_orig <- run_model_corr(simC2, analysis_omega_version="original", L=60)
# 
# ## 4) Inspecter rapidement (ex: moyennes mu)
# lapply(out_corr_deg1, function(s) apply(s$mu, 2, mean))

# ---- FIN Script 3 ----
