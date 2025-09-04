plot_clusters <- function(clusters_data) {
  
  p <- ggplot(clusters_data, aes(x = x, y = y, color = factor(cluster), fill = factor(cluster))) +
    geom_point(size = 2, alpha = 0.5) +
    scale_color_manual(values = c("1" = "#1f77b4", "2" = "#ff7f0e", "3" = "#2ca02c")) +
    scale_fill_manual(values = c("1" = "#1f77b4", "2" = "#ff7f0e", "3" = "#2ca02c")) +
    facet_grid(season ~ scenario) +
    theme_minimal() +
    labs(title = paste("Clusters"),
         x = "Coordonnée X", 
         y = "Coordonnée Y",
         color = "Cluster",
         fill = "Cluster")
  
  return(p)
}
