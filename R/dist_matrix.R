dist_matrix <- function(tab) {
  n <- nrow(tab)
  D <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      dx <- tab$x[i] - tab$x[j]
      dy <- tab$y[i] - tab$y[j]
      D[i, j] <- sqrt(dx^2 + dy^2)
    }
  }
  return(D)
}
