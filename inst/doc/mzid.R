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
  xs = list(mzid_PSM = x_toy),
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
# mzid_path <- "path/to/search_results.mzid"
# 
# x_mzid <- read_dfdr_mzid(
#   mzid_path = mzid_path,
#   rank = 1L,  # competed universe: take rank-1 SpectrumIdentificationItem
# 
#   # Choose a score CV term (priority list) and interpret its direction
#   score_accession_preference = c(
#     "MS:1002257", # MS-GF:RawScore (higher is better)
#     "MS:1001330", # Mascot:score (higher is better)
#     "MS:1001328", # SEQUEST:xcorr (higher is better)
#     "MS:1002052",
#     "MS:1002049",
#     "MS:1001331",
#     "MS:1001171",
#     "MS:1001950",
#     "MS:1002466"
#   ),
#   score_direction = "auto",  # or "higher_better"/"lower_better" if auto fails
# 
#   # TDC correction: FDR_hat = (D + add_decoy)/T
#   add_decoy = 1L,
# 
#   # Strict by default: require score for all PSMs (set <1 to allow missing)
#   min_score_coverage = 1.0,
# 
#   # Fallback decoy inference if PeptideEvidence@isDecoy is not informative
#   decoy_regex = "(^##|_REVERSED$|^REV_|^DECOY_)",
# 
#   unit = "psm",
#   scope = "global",
#   provenance = list(file = basename(mzid_path))
# )
# 
# diag <- dfdr_run_all(
#   xs = list(mzid_PSM = x_mzid),
#   alpha_main = 0.01,
#   alphas = c(1e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2, 1e-1, 2e-1),
#   eps = 0.2,
#   win_rel = 0.2,
#   truncation = "warn_drop",
#   low_conf = c(0.2, 0.5)
# )
# 
# # Export outputs
# dfdr_write_report(
#   diag,
#   out_dir = "diagFDR_mzid_out",
#   formats = c("csv", "png", "manifest", "readme", "summary")
# )
# 
# # Render a single HTML report (requires rmarkdown in Suggests)
# dfdr_render_report(diag, out_dir = "diagFDR_mzid_out")

## ----pseudo-p, eval=FALSE-----------------------------------------------------
# x_mzid$p <- score_to_pvalue(
#   score = x_mzid$score,
#   method = "decoy_ecdf",
#   is_decoy = x_mzid$is_decoy
# )
# attr(x_mzid, "meta")$p_source <- "score_to_pvalue(method='decoy_ecdf' on mzid score)"
# 
# diag_p <- dfdr_run_all(
#   xs = list(mzid_PSM = x_mzid),
#   alpha_main = 0.01
# )
# 
# dfdr_write_report(
#   diag_p,
#   out_dir = "diagFDR_mzid_out_with_p",
#   formats = c("csv","png","manifest","readme","summary")
# )
# dfdr_render_report(diag_p, out_dir = "diagFDR_mzid_out_with_p")

