check_convergence <- function(samples_corr, samples_naif,
                              dim_labels = c("x_su","y_su","x_wi","y_wi"),
                              ncol_facets = 3,
                              ncol_pi = 3) {
  suppressPackageStartupMessages({
    library(ggplot2); library(dplyr); library(tidyr); library(stringr); library(tidyselect)
  })
  
  # --- Préparation ---
  S_corr <- as.data.frame(samples_corr)
  S_naif <- as.data.frame(samples_naif)
  
  if (!".model" %in% names(S_corr)) S_corr$.model <- "corrigé Ω"
  if (!".model" %in% names(S_naif)) S_naif$.model <- "naïf"
  if (!".chain" %in% names(S_corr)) S_corr$.chain <- 1
  if (!".chain" %in% names(S_naif)) S_naif$.chain <- 1
  
  S_corr$.iter <- seq_len(nrow(S_corr))
  S_naif$.iter <- seq_len(nrow(S_naif))
  S_corr$.group <- paste(S_corr$.model, S_corr$.chain)
  S_naif$.group <- paste(S_naif$.model, S_naif$.chain)
  
  # --- Colonnes μ et π ---
  re_mu <- "^mu\\[\\s*[0-9]+\\s*,\\s*[0-9]+\\s*\\]$"
  re_pi <- "^pi\\[\\s*[0-9]+\\s*\\]$"
  
  mu_cols <- union(grep(re_mu, names(S_corr), value=TRUE), grep(re_mu, names(S_naif), value=TRUE))
  pi_cols <- union(grep(re_pi, names(S_corr), value=TRUE), grep(re_pi, names(S_naif), value=TRUE))
  if (!length(mu_cols)) mu_cols <- grep("^mu\\[", names(S_corr), value=TRUE)
  if (!length(pi_cols)) pi_cols <- grep("^pi\\[", names(S_corr), value=TRUE)
  
  # --- Table de correspondance μ ---
  mk_mu_map <- function(cols) do.call(rbind, lapply(cols, function(nm){
    m <- stringr::str_match(nm, "^mu\\[\\s*([0-9]+)\\s*,\\s*([0-9]+)\\s*\\]$")
    data.frame(param=nm, k=as.integer(m[1,2]), d=as.integer(m[1,3]),
               label=sprintf("μ[k=%d, %s]", as.integer(m[1,2]), dim_labels[as.integer(m[1,3])]))
  }))
  mu_map <- mk_mu_map(mu_cols)
  
  # --- Long μ ---
  cols_corr_mu <- intersect(mu_cols, names(S_corr))
  cols_naif_mu <- intersect(mu_cols, names(S_naif))
  long_mu <- bind_rows(
    S_corr %>%
      dplyr::select(all_of(c(cols_corr_mu, ".iter", ".model", ".chain", ".group"))) %>%
      pivot_longer(all_of(cols_corr_mu), names_to="param", values_to="value"),
    S_naif %>%
      dplyr::select(all_of(c(cols_naif_mu, ".iter", ".model", ".chain", ".group"))) %>%
      pivot_longer(all_of(cols_naif_mu), names_to="param", values_to="value")
  ) %>% left_join(mu_map, by="param")
  
  # --- Long π ---
  cols_corr_pi <- intersect(pi_cols, names(S_corr))
  cols_naif_pi <- intersect(pi_cols, names(S_naif))
  long_pi <- bind_rows(
    S_corr %>%
      dplyr::select(all_of(c(cols_corr_pi, ".iter", ".model", ".chain", ".group"))) %>%
      pivot_longer(all_of(cols_corr_pi), names_to="param", values_to="value"),
    S_naif %>%
      dplyr::select(all_of(c(cols_naif_pi, ".iter", ".model", ".chain", ".group"))) %>%
      pivot_longer(all_of(cols_naif_pi), names_to="param", values_to="value")
  ) %>%
    mutate(k = as.integer(sub("^pi\\[\\s*([0-9]+)\\s*\\]$", "\\1", param)),
           label = sprintf("π[k=%d]", k))
  
  # --- Plots ---
  p_trace_mu <- ggplot(long_mu, aes(.iter, value, color=.group)) +
    geom_line(alpha=.7, linewidth=.4) +
    facet_wrap(~ label, scales="free_y", ncol=ncol_facets) +
    labs(title="Traceplots — μ (chaînes et modèles)", x="Itération", y="Valeur", color="Chaîne") +
    theme_minimal(base_size=11)
  
  p_dens_mu <- ggplot(long_mu, aes(value, color=.group, fill=.group)) +
    geom_density(alpha=.25, linewidth=.6) +
    facet_wrap(~ label, scales="free", ncol=ncol_facets) +
    labs(title="Densités postérieures — μ", x="Valeur", y="Densité", color="Chaîne", fill="Chaîne") +
    theme_minimal(base_size=11)
  
  p_trace_pi <- ggplot(long_pi, aes(.iter, value, color=.group)) +
    geom_line(alpha=.7, linewidth=.5) +
    facet_wrap(~ label, scales="free_y", ncol=ncol_pi) +
    labs(title="Traceplots — π (chaînes et modèles)", x="Itération", y="Valeur", color="Chaîne") +
    theme_minimal(base_size=11)
  
  p_dens_pi <- ggplot(long_pi, aes(value, color=.group, fill=.group)) +
    geom_density(alpha=.25, linewidth=.7) +
    facet_wrap(~ label, scales="free", ncol=ncol_pi) +
    labs(title="Densités postérieures — π", x="Valeur", y="Densité", color="Chaîne", fill="Chaîne") +
    theme_minimal(base_size=11)
  
  list(trace_mu = p_trace_mu, dens_mu = p_dens_mu, trace_pi = p_trace_pi, dens_pi = p_dens_pi)
}

check_convergence_from_bundle <- function(bundle, ctx, omega, ncol_facets = 3) {
  bind_chains <- function(lst, model_label) {
    do.call(rbind, lapply(seq_along(lst), function(i) {
      df <- as.data.frame(lst[[i]])
      df$.chain <- i
      df$.model <- model_label
      df
    }))
  }
  samples_corr <- bind_chains(bundle$results[[ctx]][[omega]]$corr, "corrigé Ω")
  samples_naif <- bind_chains(bundle$results[[ctx]][[omega]]$naif, "naïf")
  check_convergence(samples_corr, samples_naif, ncol_facets = ncol_facets)
}


cg <- check_convergence_from_bundle(bundle, ctx="C2", omega="deg1")
cg$trace_mu
cg$dens_mu
cg$trace_pi
cg$dens_pi



