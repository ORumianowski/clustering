# ==================================================
# GMM 4D Été/Hiver — Modèle corrigé (zeros-trick stabilisé)
# Simulation + NIMBLE + 2 chaînes MCMC + Diagnostics + Plots
# ==================================================
rm(list = ls())
set.seed(123)

suppressPackageStartupMessages({
  library(MASS)
  library(Matrix)
  library(terra)
  library(nimble)
  library(mvtnorm)
  library(ggplot2)
  library(viridisLite)
  library(patchwork)
  library(coda)    # pour diagnostics de convergence
})

# ==================================================
# 1) Simulation population GMM 4D
# ==================================================
K <- 3; D <- 4; N_tot <- 4000
pi_true <- c(0.25, 0.45, 0.30)

mu <- list(
  c(-3.5,  0.0, -3.0,  0.0),
  c( 0.0, -0.5,  0.2, -1.0),
  c( 3.5,  0.5,  3.0, -0.5)
)
for(k in 1:K) mu[[k]][4] <- mu[[k]][4] - 30  # hiver plus bas

build_sigma <- function(sd, R){
  D4 <- diag(sd,4,4); S <- D4 %*% R %*% D4
  if (!isSymmetric(S) || any(eigen(S, only.values=TRUE)$values <= 0))
    S <- as.matrix(nearPD(S)$mat)
  S
}

Sigma <- list(
  build_sigma(c(0.7,2.8,2.6,1.9),
              matrix(c(1,0,0.55,-0.15,0,1,0.10,0.30,0.55,0.10,1,0.65,-0.15,0.30,0.65,1),4,4,byrow=TRUE)),
  build_sigma(c(2.6,2.4,0.8,3.0),
              matrix(c(1,-0.7,0.4,-0.35,-0.7,1,-0.1,0.45,0.4,-0.1,1,0,-0.35,0.45,0,1),4,4,byrow=TRUE)),
  build_sigma(c(0.8,3.2,3.1,0.9),
              matrix(c(1,0,0.5,0.05,0,1,-0.2,0.25,0.5,-0.2,1,0.05,0.05,0.25,0.05,1),4,4,byrow=TRUE))
)

z <- sample(1:K, N_tot, replace=TRUE, prob=pi_true)
X4_all <- t(vapply(1:N_tot, function(i) MASS::mvrnorm(1, mu[[z[i]]], Sigma[[z[i]]]), numeric(4)))
colnames(X4_all) <- c("x_su","y_su","x_wi","y_wi")
X_su_all <- X4_all[,1:2]; X_wi_all <- X4_all[,3:4]

# ==================================================
# 2) Efforts d’échantillonnage Ω(x)
# ==================================================
gauss2d <- function(x,y,mu,Sigma){
  invS <- solve(Sigma); dx <- x-mu[1]; dy <- y-mu[2]
  Q <- invS[1,1]*dx^2 + (invS[1,2]+invS[2,1])*dx*dy + invS[2,2]*dy^2
  exp(-0.5 * Q)
}

pad <- 1
xmin_all <- min(c(X_su_all[,1], X_wi_all[,1])) - pad
xmax_all <- max(c(X_su_all[,1], X_wi_all[,1])) + pad
ymin_all <- min(c(X_su_all[,2], X_wi_all[,2])) - pad
ymax_all <- max(c(X_su_all[,2], X_wi_all[,2])) + pad

build_omega <- function(center_y){
  r <- rast(nrows=80,ncols=80,xmin=xmin_all,xmax=xmax_all,ymin=ymin_all,ymax=ymax_all)
  xy <- as.data.frame(xyFromCell(r, 1:ncell(r))); colnames(xy)<-c("x","y")
  r[] <- with(xy, 0.05 + 0.9*gauss2d(x,y,c(0,center_y),diag(c(6,3))))
  r
}

r_su <- build_omega(0)
r_wi <- build_omega(-30)
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


# ==================================================
# 8) Vérification de la convergence
# ==================================================
library(coda)

mcmc_list <- mcmc.list(as.mcmc(samples1), as.mcmc(samples2))
summ <- summary(mcmc_list)
print(summ)

diag_gelman <- gelman.diag(mcmc_list, multivariate = FALSE)
print(diag_gelman)

# Traçage des chaînes pour mu et pi
par(mfrow = c(2, 3))
traceplot(mcmc_list[, grep("^mu\\[", colnames(samples1))], main = "Traceplot μ (chaînes 1 & 2)")
traceplot(mcmc_list[, grep("^pi\\[", colnames(samples1))], main = "Traceplot π (chaînes 1 & 2)")
par(mfrow = c(1, 1))

# ==================================================
# 9) Extraction des moyennes postérieures
# ==================================================
extract_muSigma <- function(samples, K, D = 4){
  mu_cols <- grep("^mu\\[", colnames(samples), value = TRUE)
  mu_mat <- matrix(0, K, D)
  for(c in mu_cols){
    idx <- as.integer(regmatches(c, gregexpr("\\d+", c))[[1]])
    mu_mat[idx[1], idx[2]] <- mean(samples[, c])
  }
  Sigma_post <- vector("list", K)
  for(k in 1:K){
    idx <- grep(paste0("Prec\\[.*,.*,", k, "\\]"), colnames(samples))
    Prec_k <- matrix(colMeans(samples[, idx, drop = FALSE]), nrow = D)
    Sigma_post[[k]] <- solve(Prec_k)
  }
  list(mu = mu_mat, Sigma = Sigma_post)
}

post <- extract_muSigma(samples1, K)

mu_su_corr <- post$mu[, 1:2]
mu_wi_corr <- post$mu[, 3:4]
Sig_su_corr <- lapply(post$Sigma, \(S) S[1:2, 1:2])
Sig_wi_corr <- lapply(post$Sigma, \(S) S[3:4, 3:4])

mu_true_mat <- do.call(rbind, mu)
mu_su_true <- mu_true_mat[, 1:2]
mu_wi_true <- mu_true_mat[, 3:4]
Sig_su_true <- lapply(Sigma, \(S) S[1:2, 1:2])
Sig_wi_true <- lapply(Sigma, \(S) S[3:4, 3:4])


# ==================================================
# 10) Visualisations
# ==================================================

library(ggplot2)
library(viridisLite)
library(patchwork)

# ------------------ Points simulés ------------------
before_df <- rbind(
  data.frame(x = X_su_all[,1], y = X_su_all[,2], saison = "été"),
  data.frame(x = X_wi_all[,1], y = X_wi_all[,2], saison = "hiver")
)

after_df <- obs_long[, c("x", "y", "saison")]

cols_season <- c("été" = "#fdae61", "hiver" = "#4575b4")

p_before <- ggplot(before_df, aes(x, y, color = saison)) +
  geom_point(size = 0.6, alpha = 0.6) +
  scale_color_manual(values = cols_season) +
  coord_equal() + theme_minimal() +
  labs(title = "Points simulés (population complète)", x = "x", y = "y", color = "Saison")

p_after <- ggplot(after_df, aes(x, y, color = saison)) +
  geom_point(size = 0.6, alpha = 0.6) +
  scale_color_manual(values = cols_season) +
  coord_equal() + theme_minimal() +
  labs(title = "Points échantillonnés (Ω été / hiver)", x = "x", y = "y", color = "Saison")

p_points <- p_before | p_after


# ------------------ Ellipses de densité ------------------
ellipse_points <- function(mu, Sigma, r, n = 200){
  ang <- seq(0, 2*pi, length.out = n)
  circle <- cbind(cos(ang), sin(ang))
  pts <- circle %*% t(chol(Sigma)) * r
  sweep(pts, 2, mu, FUN = "+")
}

make_ellipses <- function(mu_su, Sig_su, mu_wi, Sig_wi, probs = c(0.5, 0.75, 0.9)){
  rlev <- sqrt(qchisq(probs, df = 2))
  lev <- sprintf("%d%%", round(100 * probs))
  out <- list(); idx <- 1
  for(k in 1:nrow(mu_su)){
    for(i in seq_along(rlev)){
      p1 <- ellipse_points(mu_su[k,], Sig_su[[k]], rlev[i])
      p2 <- ellipse_points(mu_wi[k,], Sig_wi[[k]], rlev[i])
      out[[idx]] <- data.frame(
        x = c(p1[,1], NA, p2[,1]),
        y = c(p1[,2], NA, p2[,2]),
        cluster = factor(k),
        level = lev[i]
      )
      idx <- idx + 1
    }
  }
  do.call(rbind, out)
}

ell_true <- make_ellipses(mu_su_true, Sig_su_true, mu_wi_true, Sig_wi_true)
ell_corr <- make_ellipses(mu_su_corr, Sig_su_corr, mu_wi_corr, Sig_wi_corr)
cols <- viridisLite::rocket(K, begin = 0.15, end = 0.8, direction = -1)

p_true <- ggplot() +
  geom_path(data = ell_true, aes(x, y, group = interaction(cluster, level),
                                 color = cluster, alpha = level), linewidth = 1.1) +
  scale_alpha_manual(values = c("50%" = 0.6, "75%" = 0.3, "90%" = 0.15)) +
  scale_color_manual(values = cols) +
  labs(title = "Ellipses VRAIES (été + hiver)", x = "x", y = "y") +
  coord_equal() + theme_minimal()

p_corr <- ggplot() +
  geom_path(data = ell_corr, aes(x, y, group = interaction(cluster, level),
                                 color = cluster, alpha = level), linewidth = 1.1) +
  scale_alpha_manual(values = c("50%" = 0.6, "75%" = 0.3, "90%" = 0.15)) +
  scale_color_manual(values = cols) +
  labs(title = "Ellipses POSTÉRIEURES (corrigé)", x = "x", y = "y") +
  coord_equal() + theme_minimal()

# ------------------ Affichage global ------------------
(p_points) / (p_true | p_corr) +
  plot_annotation(title = "Diagnostic : simulation, convergence, et estimation GMM",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))
