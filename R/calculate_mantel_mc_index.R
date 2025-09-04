calculate_mantel_mc_index <- function(df_scenario) {
  ete_data <- df_scenario %>% filter(season == "breeding")
  hiver_data <- df_scenario %>% filter(season == "wintering")
  dist_ete <- dist_matrix(ete_data)
  dist_hiver <- dist_matrix(hiver_data)
  mantel_result <- mantel(dist_ete, dist_hiver, permutations = 999)
  return(data.frame(
    Mantel_r = mantel_result$statistic,
    P_value = mantel_result$signif
  ))
}



