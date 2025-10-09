# Installer si besoin :
# install.packages("ggplot2")
# install.packages("gridExtra")
# install.packages("MASS")
# install.packages("mclust")
# install.packages("clue")

library(ggplot2)
library(gridExtra)
library(MASS)
library(clue)
library(mclust)
library(tidyverse)

set.seed(123)

# -------------------------------------------------------------------
# Fonction pour créer des clusters indépendants ou corrélés
# Description de la covariance de la gaussienne selon le cas
# -------------------------------------------------------------------
generate_data <- function(n, corr = FALSE, mean_shift = 0) {
  if (corr) {
    cov_matrix <- matrix(c(1, runif(1, -0.9, 0.9),
                           runif(1, -0.9, 0.9), 1), nrow = 2)
    cluster <- mvrnorm(n, mu = c(mean_shift, mean_shift), Sigma = cov_matrix)
  } else {
    cluster <- matrix(rnorm(2 * n), ncol = 2) + mean_shift
  }
  return(cluster)
}

# -------------------------------------------------------------------
# Générer un scénario de données 4D avec 3 clusters
# -------------------------------------------------------------------
generate_scenario <- function(n_per_cluster, corr, mean_shift, unequal_sizes = FALSE) {
  if (unequal_sizes) {
    sizes <- c(25, 50, 200)
  } else {
    sizes <- rep(n_per_cluster, 3)
  }
  
  clusters <- list()
  for (i in 1:3) {
    shift <- mean_shift * (i - 2)  # décalage de -mean_shift, 0, +mean_shift
    coords12 <- generate_data(sizes[i], corr, shift)
    coords34 <- generate_data(sizes[i], corr, shift)
    clusters[[i]] <- data.frame(
      x1 = coords12[, 1],
      x2 = coords12[, 2],
      x3 = coords34[, 1],
      x4 = coords34[, 2],
      cluster = factor(paste("Cluster", i))
    )
  }
  
  df <- do.call(rbind, clusters)
  return(df)
}

# -------------------------------------------------------------------
# Fonction pour générer les deux plots pour un scénario
# -------------------------------------------------------------------
plot_scenario <- function(df, title_suffix = "") {
  p1 <- ggplot(df, aes(x = x1, y = x2, color = cluster)) +
    geom_point(size = 1.5) +
    labs(title = paste("Plan (x1,x2)", title_suffix)) +
    theme_minimal() +
    theme(legend.position = "none")
  
  p2 <- ggplot(df, aes(x = x3, y = x4, color = cluster)) +
    geom_point(size = 1.5) +
    labs(title = paste("Plan (x3,x4)", title_suffix)) +
    theme_minimal() +
    theme(legend.position = "none")
  
  return(list(p1, p2))
}

# -------------------------------------------------------------------
# Préparer les 8 scénarios uniques (2×2×2 plan factoriel)
# -------------------------------------------------------------------
corr_values <- c(FALSE, TRUE)
sizes_values <- c(FALSE, TRUE)
barycenter_values <- c(1, 3)

n_per_cluster <- 100
scenarios <- list()
scenario_id <- 1

for (corr in corr_values) {
  for (unequal_sizes in sizes_values) {
    for (mean_shift in barycenter_values) {
      df <- generate_scenario(n_per_cluster, corr, mean_shift, unequal_sizes)
      corr_label    <- ifelse(corr, "Corrélés", "Indépendants")
      size_label    <- ifelse(unequal_sizes, "Tailles ≠", "Tailles =")
      center_label  <- ifelse(mean_shift == 1, "Centres proches", "Centres éloignés")
      
      scenarios[[scenario_id]] <- list(
        id = scenario_id,
        data = df,
        corr = corr_label,
        size = size_label,
        centers = center_label
      )
      scenario_id <- scenario_id + 1
    }
  }
}

# -------------------------------------------------------------------
# Affichage des plots pour tous les scénarios
# -------------------------------------------------------------------
all_plots <- list()
for (s in scenarios) {
  scenario_title <- paste0(
    "Scénario ", s$id, ": ",
    s$corr, " | ", s$size, " | ", s$centers
  )
  scenario_plots <- plot_scenario(s$data, title_suffix = scenario_title)
  all_plots[[length(all_plots)+1]] <- scenario_plots[[1]]  # x1-x2
  all_plots[[length(all_plots)+1]] <- scenario_plots[[2]]  # x3-x4
}

# Diviser les 16 graphiques (8 scénarios × 2 projections) en 2 groupes
group1 <- all_plots[1:8]     # Scénarios 1 à 4
group2 <- all_plots[9:16]    # Scénarios 5 à 8

grid.arrange(grobs = group1, ncol = 2, nrow = 4, top = "Scénarios 1 à 4")
grid.arrange(grobs = group2, ncol = 2, nrow = 4, top = "Scénarios 5 à 8")

# -------------------------------------------------------------------
# Fonction pour calculer la proportion de points bien classés
# -------------------------------------------------------------------
classification_accuracy <- function(true_labels, predicted_labels) {
  cm <- table(true_labels, predicted_labels)
  assignment <- clue::solve_LSAP(cm, maximum = TRUE)
  acc <- sum(cm[cbind(1:nrow(cm), assignment)]) / length(true_labels)
  return(acc)
}

# -------------------------------------------------------------------
# Résultats de classification pour chaque scénario
# -------------------------------------------------------------------
results <- data.frame(
  Scenario = integer(),
  Corr = character(),
  Sizes = character(),
  Centers = character(),
  Method = character(),
  Accuracy = numeric(),
  stringsAsFactors = FALSE
)

for (s in scenarios) {
  df <- s$data
  true_labels <- df$cluster
  X <- as.matrix(df[, 1:4])
  
  # K-means
  set.seed(123)
  km <- kmeans(X, centers = 3, nstart = 20)
  acc_km <- classification_accuracy(true_labels, km$cluster)
  
  # GMM
  gmm <- Mclust(X, G = 3)
  acc_gmm <- classification_accuracy(true_labels, gmm$classification)
  
  # Stocker les résultats
  results <- rbind(
    results,
    data.frame(
      Scenario = s$id,
      Corr = s$corr,
      Sizes = s$size,
      Centers = s$centers,
      Method = "K-means",
      Accuracy = acc_km
    ),
    data.frame(
      Scenario = s$id,
      Corr = s$corr,
      Sizes = s$size,
      Centers = s$centers,
      Method = "GMM",
      Accuracy = acc_gmm
    )
  )
}




# -------------------------------------------------------------------
# Visualisation globale des performances
# -------------------------------------------------------------------
p = ggplot(results, aes(x = Method, y = Accuracy, fill = Method)) +
  geom_boxplot(alpha = 0.6, width = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.6) +
  labs(title = "Distribution de l’efficacité par méthode",
       x = "Méthode",
       y = "Proportion bien classée") +
  scale_fill_manual(values = c("K-means" = "steelblue", "GMM" = "darkorange")) +
  theme_minimal()

p

# -------------------------------------------------------------------
# Fonction pour visualiser un scénario donné (vérité, K-means, GMM)
# -------------------------------------------------------------------
plot_clustering_results <- function(df, km_labels, gmm_labels, scenario_id) {
  # Véritables labels
  p_true <- ggplot(df, aes(x = x1, y = x2, color = cluster)) +
    geom_point(size = 1.5) +
    labs(title = paste("Scénario", scenario_id, "- Vérité")) +
    theme_minimal()
  
  # K-means
  p_km <- ggplot(df, aes(x = x1, y = x2, color = factor(km_labels))) +
    geom_point(size = 1.5) +
    labs(title = paste("Scénario", scenario_id, "- K-means")) +
    theme_minimal() +
    scale_color_discrete(name = "Cluster K-means")
  
  # GMM
  p_gmm <- ggplot(df, aes(x = x1, y = x2, color = factor(gmm_labels))) +
    geom_point(size = 1.5) +
    labs(title = paste("Scénario", scenario_id, "- GMM")) +
    theme_minimal() +
    scale_color_discrete(name = "Cluster GMM")
  
  gridExtra::grid.arrange(p_true, p_km, p_gmm, ncol = 3)
}

# Exemple : affichage pour le scénario 1
# df1 <- scenarios[[1]]$data
# X1 <- as.matrix(df1[, 1:4])
# km1 <- kmeans(X1, centers = 3, nstart = 20)
# gmm1 <- Mclust(X1, G = 3)
# plot_clustering_results(df1, km1$cluster, gmm1$classification, scenario_id = 1)


# -------------------------------------------------------------------
# Tableau 1 : différence GMM - K-means
# -------------------------------------------------------------------
diff_table <- reshape(
  results[, c("Scenario", "Corr", "Sizes", "Centers", "Method", "Accuracy")],
  idvar = c("Scenario", "Corr", "Sizes", "Centers"),
  timevar = "Method",
  direction = "wide"
) |> 
  mutate(Diff_GMM_minus_Kmeans = Accuracy.GMM - `Accuracy.K-means`) |> 
  arrange(desc(Diff_GMM_minus_Kmeans))

print("Tableau 1 : Différence d’accuracy (GMM - K-means) par scénario")
print(diff_table)

# -------------------------------------------------------------------
# Tableau 2 : comparaison moyenne des méthodes
# -------------------------------------------------------------------
# comparison <- aggregate(Accuracy ~ Method, data = results, FUN = mean)
# print("Tableau 2 : Comparaison globale des méthodes")
# print(comparison)

