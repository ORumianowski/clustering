## ==================================================
## GMM 2x2D (Été / Hiver) — version COMPILABLE (static.posix)
## Zeros-trick par saison + 2 chaînes MCMC + Diagnostics + Plots
## ==================================================

rm(list = ls())
set.seed(123)

suppressPackageStartupMessages({
  library(MASS); library(Matrix); library(terra); library(nimble)
  library(ggplot2); library(viridisLite); library(patchwork); library(coda)
})

## 0) S'aligner sur un standard C++ accepté par ta toolchain
nimble::nimbleOptions(cppStd = "c++11")

## ==================================================
## 1) Simulation population 4D
## ==================================================
K <- 3L; D <- 4L; N_tot <- 4000L
pi_true <- c(0.25, 0.45, 0.30)

mu <- list(
  c(-3.5,  0.0, -3.0,  0.0),
  c( 0.0, -0.5,  0.2, -1.0),
  c( 3.5,  0.5,  3.0, -0.5)
)
for(k in 1:K) mu[[k]][4] <- mu[[k]][4] - 30

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

## ==================================================
## 2) Effort d’échantillonnage Ω(x)
## ==================================================
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
plot(r_wi)

Omega_su_all <- terra::extract(r_su, X_su_all)[,1]
Omega_wi_all <- terra::extract(r_wi, X_wi_all)[,1]

## ==================================================
## 3) Échantillonnage (on observe une saison à la fois)
## ==================================================
sel_ete   <- runif(N_tot) < Omega_su_all
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

## ==================================================
## 4) NIMBLE — 2D MVN + logZ par saison + zeros-trick
## ==================================================

## MVN 2D stable (Cholesky) utilisable dans nimbleFunction
dmvnorm2_chol <- nimbleFunction(
  run = function(x = double(1), mean = double(1), Prec = double(2),
                 log = logical(0, default = FALSE)) {
    returnType(double(0))
    xm <- x - mean
    U <- chol(Prec)
    qf <- inprod(xm, Prec %*% xm)
    ldet <- 2.0 * log(U[1,1]) + 2.0 * log(U[2,2])  # 2D plus rapide
    logdens <- 0.5*ldet - log(2*pi) - 0.5*qf
    if (log) return(logdens) else return(exp(logdens))
  }
)

## logZ par saison (intégration sur la grille) — 2D uniquement
logZ_calc_2d <- nimbleFunction(
  run=function(mu2=double(2), Prec2=double(3), pi=double(1),
               grid=double(2), omega=double(1), A=double(0),
               K=integer(0), M=integer(0)){
    returnType(double(0))
    x2 <- numeric(2)
    m2 <- numeric(2)
    P2 <- matrix(0.0, 2, 2)
    sumZ <- 0.0
    for(m in 1:M){
      x2[1] <- grid[m,1]; x2[2] <- grid[m,2]
      mix <- 0.0
      for(k in 1:K){
        m2[1] <- mu2[k,1]; m2[2] <- mu2[k,2]
        P2[1,1] <- Prec2[1,1,k]; P2[1,2] <- Prec2[1,2,k]
        P2[2,1] <- Prec2[2,1,k]; P2[2,2] <- Prec2[2,2,k]
        mix <- mix + pi[k] * dmvnorm2_chol(x2, m2, P2, FALSE)
      }
      sumZ <- sumZ + omega[m] * mix
    }
    Z <- A * (sumZ / M)
    if(Z < 1e-300) Z <- 1e-300
    return(log(Z))
  }
)

## Modèle : on évalue la vraisemblance 2D de la saison observée uniquement
code_mod <- nimbleCode({
  
  ## Normalisations
  logZ_su <- logZ_calc_2d(mu_su[1:K,1:2], Prec_su[1:2,1:2,1:K], pi[1:K],
                          grid_su[1:M_su,1:2], omega_su_grid[1:M_su], A_su, K, M_su)
  logZ_wi <- logZ_calc_2d(mu_wi[1:K,1:2], Prec_wi[1:2,1:2,1:K], pi[1:K],
                          grid_wi[1:M_wi,1:2], omega_wi_grid[1:M_wi], A_wi, K, M_wi)
  logZ[1] <- logZ_su
  logZ[2] <- logZ_wi
  
  for(i in 1:N){
    ## Mélange été
    for(k in 1:K){
      dens_su[i,k] <- pi[k] * dmvnorm2_chol(Xsu[i,1:2], mu_su[k,1:2], Prec_su[1:2,1:2,k], FALSE)
      dens_wi[i,k] <- pi[k] * dmvnorm2_chol(Xwi[i,1:2], mu_wi[k,1:2], Prec_wi[1:2,1:2,k], FALSE)
    }
    mix_su[i] <- sum(dens_su[i,1:K])
    mix_wi[i] <- sum(dens_wi[i,1:K])
    
    ## Sélectionner la saison observée via indicateur constant ind_su[i]∈{0,1}
    mixdens[i] <- ind_su[i]*mix_su[i] + (1 - ind_su[i])*mix_wi[i]
    ll[i]      <- log(Omega[i]*mixdens[i] + 1e-300) - (ind_su[i]*logZ[1] + (1 - ind_su[i])*logZ[2])
    
    ## zeros-trick stabilisé
    ll_clip[i] <- min(max(ll[i], -700), 700)
    zeros[i]   ~ dpois(phi[i])
    phi[i]     <- pmax(-ll_clip[i] + 1e-6, 1e-6)
  }
  
  ## Priors
  pi[1:K] ~ ddirch(alpha[1:K])
  
  for(k in 1:K){
    mu_su[k,1:2] ~ dmnorm(mu0_2[1:2], prec=Prec0_2[1:2,1:2])
    mu_wi[k,1:2] ~ dmnorm(mu0_2[1:2], prec=Prec0_2[1:2,1:2])
    Prec_su[1:2,1:2,k] ~ dwish(R2[1:2,1:2], df2)
    Prec_wi[1:2,1:2,k] ~ dwish(R2[1:2,1:2], df2)
  }
})

## ==================================================
## 5) Données/constants (2D par saison) + types OK
## ==================================================
X <- as.matrix(obs_long[,c("x_su","y_su","x_wi","y_wi")]); storage.mode(X) <- "double"
Xsu <- X[,1:2,drop=FALSE]; Xwi <- X[,3:4,drop=FALSE]
Omega <- as.numeric(obs_long$Omega)
ss <- ifelse(obs_long$saison=="été", 1L, 2L)
ind_su <- as.numeric(ss == 1L)   # indicateur CONSTANT (0/1)
zeros <- rep(0L,N)

su_df <- as.data.frame(r_su, xy=TRUE); names(su_df)[3] <- "omega"
wi_df <- as.data.frame(r_wi, xy=TRUE); names(wi_df)[3] <- "omega"
A_su <- as.numeric(prod(res(r_su))); A_wi <- as.numeric(prod(res(r_wi)))
M_su <- 2000L; M_wi <- 2000L  # ↑ augmenter à 4000L après premier run

set.seed(42)
grid_su <- as.matrix(su_df[sample.int(nrow(su_df), M_su, TRUE), c("x","y")])
grid_wi <- as.matrix(wi_df[sample.int(nrow(wi_df), M_wi, TRUE), c("x","y")])
storage.mode(grid_su) <- "double"; storage.mode(grid_wi) <- "double"
omega_su_grid <- as.numeric(su_df$omega[sample.int(nrow(su_df), M_su, TRUE)])
omega_wi_grid <- as.numeric(wi_df$omega[sample.int(nrow(wi_df), M_wi, TRUE)])

## Hyper-priors 2D
mu0_2   <- c(0,0)
Prec0_2 <- diag(2)*1e-2
R2      <- diag(2)
df2     <- 2L + 2L

## Inits
bar_su <- colMeans(Xsu); bar_wi <- colMeans(Xwi)
mu_su_init <- matrix(rep(bar_su, each=K), K, 2) + matrix(rnorm(K*2,0,1),K,2)
mu_wi_init <- matrix(rep(bar_wi, each=K), K, 2) + matrix(rnorm(K*2,0,1),K,2)
Prec_su_init <- array(0, dim=c(2,2,K)); for(k in 1:K) Prec_su_init[,,k] <- diag(2)*0.5
Prec_wi_init <- array(0, dim=c(2,2,K)); for(k in 1:K) Prec_wi_init[,,k] <- diag(2)*0.5
inits1 <- list(mu_su=mu_su_init, mu_wi=mu_wi_init,
               Prec_su=Prec_su_init, Prec_wi=Prec_wi_init,
               pi=rep(1/K,K))
inits2 <- list(mu_su=mu_su_init + matrix(rnorm(K*2,0,0.5),K,2),
               mu_wi=mu_wi_init + matrix(rnorm(K*2,0,0.5),K,2),
               Prec_su=Prec_su_init, Prec_wi=Prec_wi_init,
               pi=rep(1/K,K))

## Constantes & data
constants <- list(
  N=N, K=K,
  M_su=M_su, M_wi=M_wi,
  A_su=A_su, A_wi=A_wi,
  alpha=rep(1,K),
  mu0_2=mu0_2, Prec0_2=Prec0_2,
  R2=R2, df2=df2,
  ind_su = ind_su
)
data_list <- list(
  Xsu = Xsu, Xwi = Xwi,
  Omega = Omega, zeros = zeros,
  grid_su=grid_su, grid_wi=grid_wi,
  omega_su_grid=omega_su_grid, omega_wi_grid=omega_wi_grid
)

## ==================================================
## 6) Compilation & MCMC (compilé)
## ==================================================
## 1) Construire le modèle (comme tu fais déjà)
model <- nimbleModel(code_mod,
                     data = data_list,
                     constants = constants,  # mets ss dans constants !
                     inits = inits1,
                     check = FALSE)

## 2) Construire le MCMC
conf <- configureMCMC(model, monitors = c("mu","pi","Prec"), useConjugacy = FALSE)
mcmc <- buildMCMC(conf)

## 3) Lancer SANS compilation (interprété)
cat("➡️ Chaîne 1 (INTERPRÉTÉE) optimisée...\n")
samples1 <- runMCMC(mcmc, niter=2000, nburnin=500, thin=2, inits=inits1, setSeed=100)

cat("➡️ Chaîne 2 (INTERPRÉTÉE) optimisée...\n")
samples2 <- runMCMC(mcmc, niter=2000, nburnin=500, thin=2, inits=inits2, setSeed=200)

## ==================================================
## 7) Diagnostics rapides
## ==================================================
suppressPackageStartupMessages({
  library(MASS); library(Matrix); library(terra); library(nimble)
  library(ggplot2); library(viridisLite); library(patchwork); library(coda)
})

suppressPackageStartupMessages({
  library(coda)
})

# Harmoniser au cas où runMCMC renvoie data.frame vs matrix
to_mcmc <- function(smp) {
  if (is.data.frame(smp)) smp <- as.matrix(smp)
  coda::mcmc(smp)
}
mcmc_list <- coda::mcmc.list(to_mcmc(samples1), to_mcmc(samples2))

print(summary(mcmc_list))
print(gelman.diag(mcmc_list, multivariate = FALSE))

# Traceplots rapides
par(mfrow = c(2, 3))
traceplot(mcmc_list[, grep("^mu\\[", colnames(as.matrix(samples1)))], main = "Trace μ")
traceplot(mcmc_list[, grep("^pi\\[", colnames(as.matrix(samples1)))], main = "Trace π")
par(mfrow = c(1, 1))

# ==================================================
# 9) Noms lisibles + traceplots/densités par paramètre
# ==================================================
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
