## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4.5
)

## ----toy----------------------------------------------------------------------
library(diagFDR)

set.seed(10)

n <- 8000
toy <- data.frame(
  id = as.character(seq_len(n)),
  is_decoy = sample(c(FALSE, TRUE), n, replace = TRUE, prob = c(0.98, 0.02)),
  run = sample(paste0("run", 1:3), n, replace = TRUE),
  score = c(rnorm(n * 0.98, mean = 7, sd = 1), rnorm(n * 0.02, mean = 5, sd = 1))[seq_len(n)],
  pep = NA_real_
)

# Reconstruct q-values by simple TDC from score (for toy data only):
# Here we mimic a typical "higher score is better" setting.
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

diag <- dfdr_run_all(xs = list(MaxQuant_PSM = x_toy), alpha_main = 0.01)

## ----headline-----------------------------------------------------------------
diag$tables$headline

## ----stability-plots----------------------------------------------------------
diag$plots$dalpha
diag$plots$cv

## ----boundary-stability-------------------------------------------------------
diag$plots$dwin
diag$plots$elasticity

## ----equal-chance-------------------------------------------------------------
diag$tables$equal_chance_pooled
diag$plots$equal_chance__MaxQuant_PSM

## ----real-maxquant, eval=FALSE------------------------------------------------
# library(diagFDR)
# 
# mq_path <- "path/to/msms.txt"
# 
# # Read msms.txt and reconstruct TDC q-values using MaxQuant Score and Reverse indicator.
# # - Reverse == "+" is treated as a decoy indicator.
# # - Score is assumed "higher is better".
# # - q-values are computed using FDR(i) = (D(i)+add_decoy)/T(i) and q(i)=min_{j>=i} FDR(j).
# x_mq <- read_dfdr_maxquant_msms(
#   path = mq_path,
#   pep_mode = "sanitize",          # or "drop" if PEP contains values >1
#   exclude_contaminants = TRUE,
#   add_decoy = 1L,
#   unit = "psm",
#   scope = "global",
#   provenance = list(tool = "MaxQuant", file = basename(mq_path))
# )
# 
# # Run diagnostics
# diag <- dfdr_run_all(
#   xs = list(MaxQuant_PSM = x_mq),
#   alpha_main = 0.01,
#   alphas = c(1e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2, 1e-1, 2e-1),
#   eps = 0.2,
#   win_rel = 0.2,
#   truncation = "warn_drop",
#   low_conf = c(0.2, 0.5)
# )
# 
# # Inspect headline diagnostics
# diag$tables$headline
# 
# # Export tables/plots + human-readable summary
# dfdr_write_report(
#   diag,
#   out_dir = "diagFDR_maxquant_out",
#   formats = c("csv", "png", "manifest", "readme", "summary")
# )
# 
# # Optional: render a single HTML report (requires rmarkdown in Suggests)
# dfdr_render_report(diag, out_dir = "diagFDR_maxquant_out")
# 
# 
# library(diagFDR)
# 
# mq_path <- "path/to/msms.txt"
# 
# # Read msms.txt and reconstruct TDC q-values
# x_mq <- read_dfdr_maxquant_msms(
#   path = mq_path,
#   pep_mode = "sanitize",
#   exclude_contaminants = TRUE,
#   add_decoy = 1L,
#   unit = "psm",
#   scope = "global",
#   provenance = list(tool = "MaxQuant", file = basename(mq_path))
# )
# 
# # Run diagnostics with automatic pseudo p-value computation
# # This will use score_to_pvalue(method="decoy_ecdf") internally
# diag <- dfdr_run_all(
#   xs = list(MaxQuant_PSM = x_mq),
#   alpha_main = 0.01,
#   compute_pseudo_pvalues = TRUE,   # generates pseudo p-values from scores
#   pseudo_pvalue_method = "decoy_ecdf",  # Most defensible method for arbitrary scores
#   p_stratify = "run"  # Optional: stratify p-value diagnostics by run
# )
# 
# # diag will contain:
# # - All standard target-decoy diagnostics
# # - PLUS p-value calibration plots and tables
# # - PLUS Storey pi0 diagnostics
# # - PLUS BH comparison diagnostics
# 
# # Inspect headline + p-value calibration
# diag$tables$headline
# diag$tables$p_calibration_summary
# 
# # Export everything
# dfdr_write_report(
#   diag,
#   out_dir = "diagFDR_maxquant_out",
#   formats = c("csv", "png", "manifest", "readme", "summary")
# )
# 
# # Render HTML report (will include p-value diagnostics section)
# dfdr_render_report(diag, out_dir = "diagFDR_maxquant_out")

## ----pseudo-p, eval=FALSE-----------------------------------------------------
# # Create pseudo-p-values from the decoy score tail and rerun diagnostics with p-value checks.
# x_mq$p <- score_to_pvalue(
#   score = x_mq$score,
#   method = "decoy_ecdf",
#   is_decoy = x_mq$is_decoy
# )
# attr(x_mq, "meta")$p_source <- "score_to_pvalue(method='decoy_ecdf' on MaxQuant Score)"
# 
# diag_p <- dfdr_run_all(
#   xs = list(MaxQuant_PSM = x_mq),
#   alpha_main = 0.01,
#   p_stratify = "run"   # optional stratification if run column is meaningful
# )
# 
# dfdr_write_report(diag_p, out_dir = "diagFDR_maxquant_out_with_p", formats = c("csv","png","manifest","readme","summary"))
# dfdr_render_report(diag_p, out_dir = "diagFDR_maxquant_out_with_p")

