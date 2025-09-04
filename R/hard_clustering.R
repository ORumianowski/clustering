# library(cluster)
library(dplyr)
library(ggplot2)
library(tidyr)


# load(file = "donnees_connectivite.RData")



# -----------------------------
# Étape 4 : plot des clusters et sous-clusters
# -----------------------------
# On récupère les coordonnées moyenne été + hiver pour placer les points
#
# library(dplyr)
# library(dplyr)
#
# # On filtre uniquement le scénario "intermediaire"
# df_scenario <- donnees_connectivite %>%
#   filter(scenario == "intermediaire") %>%
#   arrange(saison, id)  # S'assure que l'ordre correspond à final_clusters
#
# # On ajoute les clusters estimés
# df_with_clusters <- df_scenario %>%
#   mutate(global_id = row_number()) %>%  # identifiant unique pour la fusion
#   left_join(final_clusters %>%
#               rename(global_id = id,
#                      cluster_estime = cluster,
#                      subcluster_estime = subcluster),
#             by = "global_id") %>%
#   select(-global_id)  # on supprime l'identifiant temporaire
#
# # Résultat final
# df_with_clusters
