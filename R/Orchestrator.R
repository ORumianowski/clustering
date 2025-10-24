# =============================================================
# Script 4/4 — Orchestrateur (18 runs = 3 contextes × 3 Ω d'analyse × 2 modèles)
# =============================================================
# Réutilise:
#   - Script 1: contextes + Ω
#   - Script 2: simulation (toujours Ω original)
#   - Script 3: analyse NIMBLE (naïf & corrigé) + runners 2 chaînes
# 
# Ce script:
#   • Simule 1 dataset par contexte (C1,C2,C3) avec Ω original
#   • Pour chaque contexte, lance:
#       - Modèle NAÏF (2 chaînes) — répété pour étiqueter les 3 versions d'Ω d'analyse
#       - Modèle CORRIGÉ (2 chaînes) avec Ω d'analyse ∈ {original, deg1, deg2}
#   • Stocke tout dans une liste hiérarchique: bundle$sims / bundle$results
#   • Fournit des helpers pour résumer les sorties
# =============================================================

suppressPackageStartupMessages({
  library(nimble)
})

# --- Sourcing des scripts précédents ---
if (!exists("get_context_params")) source("R/1_contexts_omegas.R")
if (!exists("simulate_dataset")) source("R/2_simulate_datasets.R")
if (!exists("run_model_naif")) source("R/3_models_analysis.R")

# ---------- Helpers de résumé ----------
# Retourne une liste avec les moyennes a posteriori de mu par chaîne
summarize_mu_means <- function(samples_list){
  lapply(samples_list, function(s)
    if (!is.null(s$mu)) {
      # s est un coda::mcmc ? runMCMC retourne data.frame; mu.* colonnes
      nm <- grep("^mu\\[", colnames(s), value = TRUE)
      if (length(nm)==0) return(NULL)
      colMeans(s[, nm, drop=FALSE])
    } else NULL)
}

# ---------- Orchestrateur principal ----------
run_all <- function(N_tot=5000, n_pix_side=60,
                    niter=3000, nburnin=1000, thin=5,
                    seed_sim_base=100, L_corr=60,
                    contexts=c("C1","C2","C3"),
                    omega_analysis=c("original","deg1","deg2")){
  
  # 1) SIMULATIONS (Ω ORIGINAL) — 1 dataset par contexte
  sims <- list()
  for (i in seq_along(contexts)){
    ctx <- contexts[i]
    message(sprintf("[SIM] Contexte %s (Ω simulation: original)", ctx))
    sims[[ctx]] <- simulate_dataset(ctx, N_tot=N_tot, seed=seed_sim_base + i, n_pix_side=n_pix_side)
  }
  
  # 2) ANALYSES (18 runs -> optimisé: naïf UNE FOIS par contexte)
  results <- list()
  for (ctx in contexts){
    results[[ctx]] <- list()
    sim <- sims[[ctx]]
    
    # --- NAÏF: une seule exécution (2 chaînes) par contexte ---
    message(sprintf("[RUN] %s — modèle naïf (2 chaînes, exécuté une seule fois)", ctx))
    res_naif_once <- run_model_naif(sim, niter=niter, nburnin=nburnin, thin=thin)
    
    # --- CORRIGÉ: une exécution par Ω d'analyse ---
    for (om in omega_analysis){
      key <- paste(ctx, om, sep="__")
      message(sprintf("[RUN] %s — modèle corrigé Ω=%s (2 chaînes)", ctx, om))
      res_corr <- run_model_corr(sim, analysis_omega_version = om, L=L_corr,
                                 niter=niter, nburnin=nburnin, thin=thin)
      # On conserve la compatibilité de structure: ranger le même résultat naïf avec chaque entrée Ω
      results[[ctx]][[om]] <- list(naif = res_naif_once, corr = res_corr)
    }
  }
  
  list(sims=sims, results=results, meta=list(N_tot=N_tot, niter=niter, nburnin=nburnin, thin=thin, L_corr=L_corr,
                                             contexts=contexts, omega_analysis=omega_analysis,
                                             note="Naïf exécuté une seule fois par contexte et réutilisé pour chaque Ω"))
}

# ---------- Exemples d'utilisation (décommenter si besoin) ----------
bundle <- run_all(N_tot=3000, n_pix_side=60, niter=200, nburnin=1, thin=1, seed_sim_base=100, L_corr=60)

# Accéder aux résultats: ex. moyennes de mu pour C2 / Ω deg1 / corrigé
mu_means_chains <- summarize_mu_means(bundle$results[["C2"]][["deg1"]]$corr)
str(mu_means_chains)

# Sauvegarde
saveRDS(bundle, file = sprintf("gmm4d_omega_18runs_%s.rds", format(Sys.time(), "%Y%m%d_%H%M%S")))

# ---- FIN Script 4 ----
