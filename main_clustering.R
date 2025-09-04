library(cluster)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(tidyr)
library(vegan)     
library(ggforce) 

devtools::load_all("R")



# Simulate clusters -------------------------------------------------------

# Set parameters
nDots <- 30
nClusters <- 3

# Generate breeding season clusters
breeding_clusters <- generate_clusters_breeding(n_points = nDots) |> 
  mutate(season = "breeding")

# Define scenarios
scenarios <- c("forte", "intermediaire", "faible")

# Generate wintering clusters for each scenario
wintering_clusters <- lapply(scenarios, function(scenario) {
  generate_clusters_3types(n_points = 30, scenario) |> 
    mutate(season = "wintering")
}) %>% 
  bind_rows()

# Duplicate breeding clusters for each scenario
breeding_clusters <- lapply(scenarios, function(scenario) {
  breeding_clusters %>% 
    mutate(scenario = scenario)
}) %>% 
  bind_rows()

# Combine breeding and wintering data
clusters <- rbind(breeding_clusters, wintering_clusters)

# Convert scenario to factor with specified order to maintain the order
clusters$scenario <- factor(clusters$scenario, levels = scenarios)

plot_clusters(clusters)

# Mantel index ------------------------------------------------------------

calculate_mantel_mc_index(clusters|> 
                            subset(scenario == "faible"))
calculate_mantel_mc_index(clusters|> 
                            subset(scenario == "intermediaire"))
calculate_mantel_mc_index(clusters|> 
                            subset(scenario == "forte"))




plots <- lapply(scenarios, function(scenario) {
  create_connectivity_plot(scenario, clusters)
})

final_plot <- grid.arrange(grobs = plots, ncol = 3)
print(final_plot)


# soft --------------------------------------------------------------------



clusters = generate_clusters_2D(n_points = 30, n_clusters = 3) |> 
  select(x, y)




# --- Gaussian density ---
# --- Gaussian density ---
dmvnorm <- function(x, mean, cov) {
  k <- length(mean)
  det_cov <- det(cov)
  if (det_cov <= 0) {
    cov <- cov + diag(1e-6, k) # régularisation
    det_cov <- det(cov)
  }
  inv_cov <- solve(cov)
  norm_const <- 1 / ((2 * pi)^(k/2) * sqrt(det_cov))
  exp_term <- exp(-0.5 * rowSums((x - mean) %*% inv_cov * (x - mean)))
  return(norm_const * exp_term)
}

# --- Expectation-Maximization pour GMM ---
gmm_em <- function(data, K = 3, max_iter = 100, tol = 1e-6) {
  N <- nrow(data)
  d <- ncol(data)
  
  # Initialisation aléatoire
  set.seed(123)
  means <- data[sample(1:N, K), ]
  covariances <- array(0, c(d, d, K))
  for (k in 1:K) covariances[,,k] <- diag(d)
  weights <- rep(1/K, K)
  
  log_likelihood_old <- -Inf
  
  for (iter in 1:max_iter) {
    # --- E-step ---
    resp <- matrix(0, N, K)
    for (k in 1:K) {
      resp[,k] <- weights[k] * dmvnorm(data, means[k,], covariances[,,k])
    }
    resp <- resp / rowSums(resp)
    
    # --- M-step ---
    Nk <- colSums(resp)
    for (k in 1:K) {
      if (Nk[k] < 1e-6) next  # éviter une composante vide
      means[k,] <- colSums(resp[,k] * data) / Nk[k]
      diff <- data - matrix(means[k,], nrow=N, ncol=d, byrow=TRUE)
      covariances[,,k] <- t(diff) %*% (diff * resp[,k]) / Nk[k] + diag(1e-6, d)
      weights[k] <- Nk[k] / N
    }
    
    # --- Log-likelihood ---
    log_likelihood <- sum(log(rowSums(sapply(1:K, function(k) {
      weights[k] * dmvnorm(data, means[k,], covariances[,,k])
    }))))
    
    if (abs(log_likelihood - log_likelihood_old) < tol) {
      break
    }
    log_likelihood_old <- log_likelihood
  }
  
  list(means = means, covariances = covariances, weights = weights,
       responsibilities = resp, logLik = log_likelihood, iter = iter)
}

# Appliquer GMM explicite
res <- gmm_em(as.matrix(clusters), K = 3)

# Cluster final = argmax des responsabilités
cl <- apply(res$responsibilities, 1, which.max)

# Visualisation
plot(clusters$x, clusters$y)
plot(clusters$x, clusters$y, col = cl, pch = 19)


