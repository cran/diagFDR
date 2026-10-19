# diagFDR 0.1.2

## Diagnostics

* `dfdr_ipe()`: decoy-reference PEP concordance (internal PEP calibration
  error, IPE) with the target-decoy comparator `D_b/T_b` (mean target PEP per
  bin, weights `n_b = T_b + D_b`, bins with mean target PEP <= 0.5 and
  `n_b >= 200`). The headline `IPE` uses it; `IPE_pooled` gives the same
  summary with the pooled comparator `D_b/(T_b + D_b)`.
* `dfdr_headline()` returns the decoy-support index `support_index`
  (`D^-1/2`) and the FDP-CV benchmark `CV_benchmark` (`sqrt(2) * D^-1/2`).
* `flag_headline()`: review alerts, not a verdict. Flag values `OK`,
  `CAUTION*`, `ALERT_*`, `NE_*` (not evaluable) and `NA`; status `NO_ALERT`,
  `CAUTION` or `ALERT`. IPE alerts above 0.05. The equal-chance flag is
  `NE_two_pass` on two-pass exports and `CAUTION_imbalance_pass_unknown` on a
  DIA-NN list whose pass status was not given.
* `dfdr_alert_basis()`: evidential basis of each alert (derived check,
  simulation-informed alert, heuristic alert, descriptive).
* `dfdr_p_summary()` / `dfdr_p_interpretation()`: pi0 by the eight cp4p
  estimators on targets + decoys, targets and decoys, grouped into a converging
  set (consensus value) and the others; p-value threshold, and score or PEP
  threshold, giving an estimated FDR of alpha on targets with the consensus
  pi0 and with pi0 = 1 (Benjamini-Hochberg); comparison with target-decoy
  counting; pooled targets + decoys check; caveats.
* `dfdr_pep_summary()` / `dfdr_pep_interpretation()`: whether PEPs agree with,
  under- or overestimate the decoy-derived reference, and the sum of the PEPs
  of the accepted targets against the number of accepted decoys.

## Inputs

* `as_dfdr_tbl()` and the DIA-NN builders take `two_pass` (set `TRUE` for
  `--reanalyse` / match-between-runs reports); the DIA-NN builders say once
  per session when it is not given.
* `dfdr_run_all(compute_pseudo_pvalues = TRUE)` warns, and records the reason
  in the list metadata (`p_source`), when a list has no `score` column (DIA-NN
  2.x reports carry no `CScore`). `pseudo_pvalue_score = "pep"` uses `-PEP` as
  an approximate surrogate score.
* `mokapot_competed_universe(q_mode = c("as_exported", "recount"))`: mokapot
  computes q-values over all PSMs it receives; `"recount"` recomputes
  `(D+1)/T` q-values on the competed winners.
* DIA-NN builders keep the largest `score` of a collapsed unit (higher is
  better), and the smallest q-value and PEP.

## HTML report

* Executive summary with alerts to review, each tagged with its evidential
  basis; lists without alert shown as "no internal alert"; checks that are
  not evaluable shown as "n/e".
* Every table is followed by a legend of its columns; small values are shown
  in scientific notation; formulas are written in HTML (no MathJax needed).
* "What the p-value diagnostics say" and a PEP interpretation, written from
  the summaries above.
* The recomputed decoy ratio is labelled `D/T`; the Jaccard guide treats
  overlap as descriptive.
* `dfdr_render_report()` accepts a relative `out_dir` and renders from a
  temporary copy of the template.
* A plot that cannot be built for a list (e.g. too few p-values) is skipped
  with a warning.
