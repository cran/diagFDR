## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4.5
)

## ----toy-data-----------------------------------------------------------------
library(diagFDR)

set.seed(1)

n <- 3000
toy_global <- data.frame(
  id = paste0("P", seq_len(n)),
  is_decoy = sample(c(FALSE, TRUE), n, replace = TRUE, prob = c(0.97, 0.03)),
  q = pmin(1, runif(n)^3),       # skew toward small q-values
  pep = NA_real_,
  run = NA_character_,
  score = NA_real_
)

x_global <- as_dfdr_tbl(
  toy_global,
  unit = "precursor",
  scope = "global",
  q_source = "toy",
  q_max_export = 0.5
)

diag <- dfdr_run_all(
  xs = list(global = x_global),
  alpha_main = 0.01,
  alphas = c(1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2, 1e-1, 2e-1),
  low_conf = c(0.2, 0.5)
)

## ----headline-----------------------------------------------------------------
diag$tables$headline

## ----plots-stability----------------------------------------------------------
diag$plots$dalpha
diag$plots$cv

## ----plot-dwin----------------------------------------------------------------
diag$plots$dwin

## ----plot-elasticity----------------------------------------------------------
diag$plots$elasticity

## ----equal-chance-------------------------------------------------------------
diag$tables$equal_chance_pooled
diag$plots$equal_chance__global

## ----real-diann, eval=FALSE---------------------------------------------------
# # Requires arrow
# rep <- read_diann_parquet("path/to/report.parquet")
# 
# # (A) Global precursor list using Global.Q.Value
# # Recommended for experiment-wide (pooled) lists.
# x_global_gq <- diann_global_precursor(
#   rep,
#   q_col = "Global.Q.Value",
#   q_max_export = 0.5,
#   unit = "precursor",
#   scope = "global",
#   q_source = "Global.Q.Value"
# )
# 
# # (B) Run×precursor universe using run-wise Q.Value
# # Recommended for per-run decisions / QC.
# x_runx <- diann_runxprecursor(
#   rep,
#   q_col = "Q.Value",
#   q_max_export = 0.5,
#   id_mode = "runxid",
#   unit = "runxprecursor",
#   scope = "runwise",
#   q_source = "Q.Value"
# )
# 
# # (C) Scope misuse comparator: min run-wise q over runs per precursor (anti-pattern)
# # Useful for demonstrating/diagnosing scope mismatch.
# x_minrun <- diann_global_minrunq(
#   rep,
#   q_col = "Q.Value",
#   q_max_export = 0.5,
#   unit = "precursor",
#   scope = "aggregated",
#   q_source = "min_run(Q.Value)"
# )
# 
# diag <- dfdr_run_all(
#   xs = list(global = x_global_gq, runx = x_runx, minrun = x_minrun),
#   alpha_main = 0.01,
#   compute_pseudo_pvalues = TRUE  # <-- This adds p-value diagnostics
# )
# 
# # Compare accepted lists across scopes (Jaccard overlap across alpha)
# scope_tbl <- dfdr_scope_disagreement(
#   x1 = x_global_gq,
#   x2 = x_minrun,
#   alphas = c(1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2),
#   label1 = "Global.Q.Value",
#   label2 = "min_run(Q.Value)"
# )
# 
# # Write outputs to disk (tables + plots; optionally PPTX)
# dfdr_write_report(diag, out_dir = "diagFDR_diann_out", formats = c("csv", "png", "manifest", "readme", "summary"))
# 
# # Render a single HTML report (requires rmarkdown in Suggests)
# dfdr_render_report(diag, out_dir = "diagFDR_diann_out")

