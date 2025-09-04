

create_connectivity_plot <- function(scenario, clusters_data) {
  
  couleurs <- c("darkred", "darkorange", "darkgreen")
  
  # Filtrer les données pour le scénario
  df <- clusters_data %>% filter(scenario == !!scenario)
  
  # Préparer segments été ↔ hiver
  liens <- df %>%
    tidyr::pivot_wider(names_from = season, values_from = c(x, y),
                       id_cols = c(cluster, id))
  
  # RECALCULER l'indice de Mantel pour ce scénario spécifique
  mantel_result <- calculate_mantel_mc_index(df)
  
  # Calcul des ellipses (moyenne, cov et niveaux de confiance)
  ellipses <- df %>%
    group_by(cluster, season) %>%
    do({
      mu <- colMeans(cbind(.$x, .$y))
      covm <- cov(cbind(.$x, .$y))
      data.frame(
        x0 = mu[1],
        y0 = mu[2],
        a50 = sqrt(qchisq(0.5, df = 2)),
        a90 = sqrt(qchisq(0.9, df = 2)),
        covm = I(list(covm))
      )
    }) %>% ungroup()
  
  # Plot
  p <- ggplot(df, aes(color = factor(cluster), fill = factor(cluster))) +
    geom_segment(data = liens,
                 aes(x = x_breeding, y = y_breeding, xend = x_wintering, yend = y_wintering,
                     color = factor(cluster)),
                 linewidth = 0.6, alpha = 0.4, inherit.aes = FALSE) +
    # Ellipses 50%
    geom_ellipse(data = ellipses,
                 aes(x0 = x0, y0 = y0, a = a50, b = a50, angle = 0,
                     fill = factor(cluster)),
                 inherit.aes = FALSE, alpha = 0.2, color = NA) +
    # Ellipses 90%
    geom_ellipse(data = ellipses,
                 aes(x0 = x0, y0 = y0, a = a90, b = a90, angle = 0,
                     fill = factor(cluster)),
                 inherit.aes = FALSE, alpha = 0.1, color = NA) +
    geom_point(aes(x = x, y = y), size = 2, alpha = 0.7) +
    scale_color_manual(values = couleurs, name = "Cluster") +
    scale_fill_manual(values = couleurs, name = "Cluster") +
    theme_minimal(base_size = 14) +
    labs(title = paste("Connectivité", scenario),
         x = "Longitude", y = "Latitude") +
    geom_text(
      aes(x = min(df$x) + (max(df$x) - min(df$x)) * 0.1, 
          y = max(df$y) - (max(df$y) - min(df$y)) * 0.1, 
          label = paste("Mantel r =", round(mantel_result$Mantel_r, 3),
                        "\np-value =", round(mantel_result$P_value, 4))),
      color = "black", size = 4, vjust = 1, hjust = 0, inherit.aes = FALSE
    )
  
  return(p)
}

