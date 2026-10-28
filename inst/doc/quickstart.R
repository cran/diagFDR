## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7, fig.height = 4)

## -----------------------------------------------------------------------------
library(diagFDR)
sim <- sim_tdc_example          # 20,000 units, Gaussian nulls, correct specification
sim

## -----------------------------------------------------------------------------
h <- dfdr_headline(sim, alpha = 0.01)
h[, c("alpha", "T_alpha", "D_alpha", "FDR_hat", "support_index", "CV_benchmark")]

## -----------------------------------------------------------------------------
dfdr_realised_fdp(sim, alpha = c(0.05, 0.01, 0.002))

## -----------------------------------------------------------------------------
dfdr_check_cv_benchmark(B = 30, n = 20000, alpha = 0.01)[, c("alpha", "CV_FDP", "mean_D", "kappa")]

## -----------------------------------------------------------------------------
alphas <- c(0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05)
dfdr_curve_stability(sim, alphas = alphas)[, c("alpha", "T_alpha", "D_alpha", "support_index")]
dfdr_local_tail(sim, alphas = alphas, win_rel = 0.2)[, c("alpha", "D_alpha_win")]
dfdr_elasticity(sim, alphas = alphas, eps = 0.2)[, c("alpha", "jaccard")]

## -----------------------------------------------------------------------------
dfdr_ipe(sim)$IPE
dfdr_equal_chance_qbands(sim)$pooled

## -----------------------------------------------------------------------------
bad <- dfdr_simulate_tdc(n = 20000, decoy_scale = 0.85, seed = 2)
dfdr_realised_fdp(bad, alpha = 0.01)

## -----------------------------------------------------------------------------
set.seed(1)
rep <- data.frame(Protein.Group = rep(paste0("P", 1:2000), each = 3),
                  Decoy = rbinom(6000, 1, 0.12),
                  Global.PG.Q.Value = c(runif(6000, 0, 0.3)))
label_is_intrinsic(rep, "Protein.Group")
dfdr_label_rule_audit(rep, "Protein.Group", "Global.PG.Q.Value", alpha = 0.05)

## -----------------------------------------------------------------------------
pg <- diann_global_proteingroup(rep, rule = "split", two_pass = FALSE)
pk <- dfdr_picked_protein(pg, alpha = 0.05)
c(T = pk$T_alpha, D = pk$D_alpha, ratio = pk$decoy_ratio, pairs = pk$n_pairs)

## -----------------------------------------------------------------------------
truth <- attr(sim, "truth")
ent_ids <- sim$id[!sim$is_decoy & truth$is_null & runif(nrow(sim)) < 0.3]
dfdr_entrapment(sim, ent_ids, alpha = c(0.05, 0.01), r = 0.43)

## -----------------------------------------------------------------------------
a1 <- dfdr_accepted(sim, 0.01)
a2 <- dfdr_accepted(sim, 0.03)      # a looser list standing in for a mis-scoped one
d  <- dfdr_entrapment_discordant(a1, a2, ent_ids)
d
attr(d, "ratio"); attr(d, "jaccard")

