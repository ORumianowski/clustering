

generate_clusters_3types <- function(n_points = 30, scenario) {
  
  # Définition des centres selon le scénario
  if (scenario == "forte") {
    centres <- data.frame(x = c(3, 6, 9), y = c(0, 1, 2))
  } else if (scenario == "intermediaire") {
    centres <- data.frame(x = c(3, 3.8, 9), y = c(0, 0.2, 2))
  } else { # faible
    centres <- data.frame(x = c(6, 6.2, 6.1), y = c(1, 1.1, 0.9))
  }
  
  # Génération des données
  donnees_list <- lapply(1:3, function(i) {
    data.frame(
      x = rnorm(n_points, centres$x[i], 0.5),
      y = rnorm(n_points, centres$y[i], 0.5),
      cluster = i,
      id = 1:n_points,
      scenario = scenario
    )
  })
  
  return(do.call(rbind, donnees_list))
}


