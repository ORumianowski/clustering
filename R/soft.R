# Installer si nécessaire
# install.packages(c("mclust", "ggplot2", "ggforce"))

library(mclust)
library(ggplot2)
library(ggforce)

set.seed(42)
n <- 500

# --- Génération des données réelles ---
X1 <- matrix(rnorm(n*4, mean=0), ncol=4)
X2 <- matrix(rnorm(n*4, mean=0.2), ncol=4)
X3 <- matrix(rnorm(n*4, mean=10), ncol=4)

X <- rbind(X1, X2, X3)
true_labels <- rep(1:3, each=n)

# --- Appliquer GMM ---
gmm <- Mclust(X, G=3)
pred_labels <- gmm$classification


a = tibble(true_labels,
           pred_labels)

table(a)

