

generate_clusters_2D <- function(n_points = 30, n_clusters = 3) {
  centres_x <- runif(n_clusters, 0, 10)
  centres_y <- runif(n_clusters, 0, 10)
  ecart_type <- 0.5
  
  donnees_list <- lapply(1:n_clusters, function(i) {
    data.frame(
      x = rnorm(n_points, centres_x[i], ecart_type),
      y = rnorm(n_points, centres_y[i], ecart_type),
      cluster = i,
      id = 1:n_points
    )
  })
  
  return(do.call(rbind, donnees_list))
}


generate_clusters_breeding <- function(n_points = 30) {


    dt = data.frame(
      x = c(rnorm(n_points, 2, 0.5), rnorm(n_points, 5, 0.5), rnorm(n_points, 8, 0.5)),
      y = c(rnorm(n_points, 8, 0.5), rnorm(n_points, 7, 0.5), rnorm(n_points, 9, 0.5)),
      cluster = rep(1:3, each = n_points),
      id = 1:n_points
    )
  
  return(dt)
}

