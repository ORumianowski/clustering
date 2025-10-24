# ===== Bloc 3 : Convergence — traceplots & densités pour μ et π =====
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
})

check_convergence <- function(samples_corr4, samples_naif4,
                              dim_labels = c("x_su","y_su","x_wi","y_wi"),
                              ncol_facets = 3) {
  
  S_corr <- as.data.frame(samples_corr4); S_corr$.iter <- seq_len(nrow(S_corr)); S_corr$.model <- "corrigé Ω"
  S_naif <- as.data.frame(samples_naif4); S_naif$.iter <- seq_len(nrow(S_naif)); S_naif$.model <- "naïf"
  
  mu_cols <- grep("^mu\\[\\s*\\d+\\s*,\\s*\\d+\\s*\\]$", names(S_corr), value = TRUE)
  pi_cols <- grep("^pi\\[\\s*\\d+\\s*\\]$", names(S_corr), value = TRUE)
  
  mk_mu_map <- function(cols) do.call(rbind, lapply(cols, function(nm){
    m <- stringr::str_match(nm, "^mu\\[\\s*(\\d+)\\s*,\\s*(\\d+)\\s*\\]$")
    k <- as.integer(m[1,2]); d <- as.integer(m[1,3])
    data.frame(param=nm, k=k, d=d, label=sprintf("mu[k=%d, %s]", k, dim_labels[d]))
  }))
  
  mu_map <- mk_mu_map(mu_cols)
  
  long_mu <- bind_rows(
    S_corr |> select(all_of(mu_cols), .iter, .model) |> pivot_longer(all_of(mu_cols), names_to="param", values_to="value"),
    S_naif |> select(all_of(mu_cols), .iter, .model) |> pivot_longer(all_of(mu_cols), names_to="param", values_to="value")
  ) |> left_join(mu_map, by="param")
  
  long_pi <- bind_rows(
    S_corr |> select(all_of(pi_cols), .iter, .model) |> pivot_longer(all_of(pi_cols), names_to="param", values_to="value"),
    S_naif |> select(all_of(pi_cols), .iter, .model) |> pivot_longer(all_of(pi_cols), names_to="param", values_to="value")
  ) |> mutate(k = as.integer(sub("^pi\\[(\\d+)\\]$", "\\1", param)),
              label = sprintf("pi[k=%d]", k))
  
  p_trace_mu <- ggplot(long_mu, aes(.iter, value, color=.model)) +
    geom_line(alpha=.85, linewidth=.5) +
    facet_wrap(~ label, scales="free_y", ncol=ncol_facets) +
    labs(title="Traceplots — μ (par modèle)", x="itération", y="valeur", color="modèle") +
    theme_minimal(base_size=12)
  
  p_dens_mu <- ggplot(long_mu, aes(value, color=.model, fill=.model)) +
    geom_density(alpha=.25, linewidth=.7) +
    facet_wrap(~ label, scales="free", ncol=ncol_facets) +
    labs(title="Densités postérieures — μ", x="valeur", y="densité", color="modèle", fill="modèle") +
    theme_minimal(base_size=12)
  
  p_trace_pi <- ggplot(long_pi, aes(.iter, value, color=.model)) +
    geom_line(alpha=.85, linewidth=.6) +
    facet_wrap(~ label, scales="free_y") +
    labs(title="Traceplots — π", x="itération", y="valeur", color="modèle") +
    theme_minimal(base_size=12)
  
  p_dens_pi <- ggplot(long_pi, aes(value, color=.model, fill=.model)) +
    geom_density(alpha=.25, linewidth=.8) +
    facet_wrap(~ label, scales="free") +
    labs(title="Densités postérieures — π", x="valeur", y="densité", color="modèle", fill="modèle") +
    theme_minimal(base_size=12)
  
  list(trace_mu = p_trace_mu, dens_mu = p_dens_mu, trace_pi = p_trace_pi, dens_pi = p_dens_pi)
}


cg <- check_convergence(samples_corr4, samples_naif4, ncol_facets = 3)
cg$trace_mu; cg$dens_mu; cg$trace_pi; cg$dens_pi