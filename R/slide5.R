# Packages
library(ggplot2)
library(mclust)
library(ggforce)
library(gridExtra)

set.seed(42)

# --- Données simulées (groupes proches pour overlap) ---
n <- 150
x1 <- cbind(rnorm(n, -1, 0.7), rnorm(n,  0, 0.7))
x2 <- cbind(rnorm(n,  1, 0.7), rnorm(n,  0, 0.7))
x3 <- cbind(rnorm(n,  0, 0.7), rnorm(n,  1.2, 0.7))
X  <- rbind(x1, x2, x3)
df <- data.frame(x = X[,1], y = X[,2])

# --- HARD CLUSTERING (KMeans) ---
km <- kmeans(X, centers = 3, nstart = 25)
df$hard <- factor(km$cluster)

# Décalage vertical ajusté
y_shift <- 0.25

# Palette plus esthétique
my_colors <- c("#984ea3", "#4daf4a", "#ff7f00")

# Plot Hard
p1 <- ggplot(df, aes(x, y, color = hard)) +
  geom_point(size=1.8, alpha=0.9) +
  annotate("segment", x=0, y=y_shift, xend=0,    yend=-2.5+y_shift, size=1.2) +
  annotate("segment", x=0, y=y_shift, xend=-2.2, yend= 1.2+y_shift, size=1.2) +
  annotate("segment", x=0, y=y_shift, xend= 2.2, yend= 1.2+y_shift, size=1.2) +
  scale_color_manual(values=my_colors) +
  labs(title="Hard clustering (K-Means)") +
  theme_minimal() +
  theme(legend.position="none")

# --- SOFT CLUSTERING (Gaussian Mixture) ---
gmm <- Mclust(X, G=3)
df$soft <- factor(gmm$classification)

# Plot Soft
p2 <- ggplot(df, aes(x, y, color = soft)) +
  geom_point(size=1.8, alpha=0.9) +
  geom_mark_ellipse(aes(fill=soft, label=NULL), alpha=0.2, show.legend=FALSE) +
  scale_color_manual(values=my_colors) +
  scale_fill_manual(values=my_colors) +
  labs(title="Soft clustering (Gaussian Mixture)") +
  theme_minimal() +
  theme(legend.position="none")

# --- Figure finale ---
grid.arrange(p1, p2, ncol=2)
