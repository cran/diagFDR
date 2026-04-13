## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4.5
)

## ----toy----------------------------------------------------------------------
library(diagFDR)

set.seed(3)

n <- 7000
pi_decoy <- 0.03

# Decoy indicator
is_decoy <- sample(c(FALSE, TRUE), n, replace = TRUE, prob = c(1 - pi_decoy, pi_decoy))

# Targets are a mixture: some null-like (close to decoys), some true (higher score)
# This yields realistic separation and non-trivial discoveries at 1% FDR.
is_true_target <- (!is_decoy) & (runif(n) < 0.30)  # 30% of targets are "true"
is_null_target <- (!is_decoy) & (!is_true_target)

score <- numeric(n)
score[is_decoy]       <- rnorm(sum(is_decoy), mean = 0.0, sd = 1.0)
score[is_null_target] <- rnorm(sum(is_null_target), mean = 0.2, sd = 1.0)
score[is_true_target] <- rnorm(sum(is_true_target), mean = 3.0, sd = 1.0)

toy <- data.frame(
  id = paste0("psm", seq_len(n)),
  is_decoy = is_decoy,
  run = sample(paste0("run", 1:4), n, replace = TRUE),
  score = score,
  pep = NA_real_
)

# Score-based TDC q-values (higher score is better)
toy <- toy[order(toy$score, decreasing = TRUE), ]
toy$D_cum <- cumsum(toy$is_decoy)
toy$T_cum <- cumsum(!toy$is_decoy)
toy$FDR_hat <- (toy$D_cum + 1) / pmax(toy$T_cum, 1)
toy$q <- rev(cummin(rev(toy$FDR_hat)))
toy <- toy[, c("id","is_decoy","q","pep","run","score")]

x_toy <- as_dfdr_tbl(
  toy,
  unit = "psm",
  scope = "global",
  q_source = "toy TDC from score",
  q_max_export = 1,
  provenance = list(tool = "toy")
)

diag <- dfdr_run_all(
  xs = list(univ = x_toy),
  alpha_main = 0.01,
  alphas = c(1e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2, 1e-1, 2e-1),
  eps = 0.2,
  win_rel = 0.2,
  truncation = "warn_drop",
  low_conf = c(0.2, 0.5)
)

## ----headline-----------------------------------------------------------------
diag$tables$headline

if (nrow(diag$tables$headline) > 0 && diag$tables$headline$T_alpha[1] == 0) {
  cat("\nNote: No discoveries at alpha_main for this toy run. ",
      "For demonstration, consider using alpha_main = 0.02.\n", sep = "")
}

## ----stability-plots----------------------------------------------------------
diag$plots$dalpha
diag$plots$cv

## ----boundary-stability-------------------------------------------------------
diag$plots$dwin
diag$plots$elasticity

## ----equal-chance-------------------------------------------------------------
diag$tables$equal_chance_pooled
diag$plots$equal_chance__mzid_PSM

## ----real-mzid, eval=FALSE----------------------------------------------------
# library(diagFDR)
# 
# # Read efficiently (when facing big datasets)
# rep <- read_spectronaut_efficient("path/to/search_results.Report-Peptide normal.tsv",
#                                   minimal = TRUE, dec = ",")
# 
# univ_runwise <- spectronaut_runxprecursor(
#   rep,
#   q_col = "EG.Qvalue",
#   score_col = "EG.Cscore"
# )
# 
# # Run diagnostics
# diag <- dfdr_run_all(
#       list(runwise = univ_runwise),
#       compute_pseudo_pvalues=TRUE
#  )
# 
# # Export outputs
# dfdr_write_report(
#   diag,
#   out_dir = "diagFDR_spectronaut_out",
#   formats = c("csv", "png", "manifest", "readme", "summary")
# )
# 
# # Optional: render a single HTML report
# dfdr_render_report(diag, out_dir = "diagFDR_spectronaut_out")

