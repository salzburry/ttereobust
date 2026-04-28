## summarise_aiptw.R
## Performance aggregator for the AIPTW simulation results.
##
## Implements the four protocol metrics from Apr 28 methods document
## Section 6.5:
##
##   - Relative bias (%):    100 * (mean(estimate) - truth) / truth
##   - Coverage:             P(CI_lo <= truth <= CI_hi)
##   - Relative SE error:    mean(model SE) / sd(estimates) - 1
##   - Power:                P(0 not in CI)   (only meaningful when truth != 0)
##
## Reads everything under results/raw/ and results/truth/, joins them, and
## writes one tidy summary CSV per (scenario, t, target).
##
## CLI usage (from apr update/code/):
##   Rscript summarise_aiptw.R
##   Rscript summarise_aiptw.R --out my_summary.csv

suppressPackageStartupMessages({
  library(dplyr)
})

CODE_DIR    <- getwd()
RESULTS_DIR <- file.path(CODE_DIR, "results")
RAW_DIR     <- file.path(RESULTS_DIR, "raw")
TRUTH_DIR   <- file.path(RESULTS_DIR, "truth")


parse_cli_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list(out = file.path(RESULTS_DIR, "summary.csv"))
  i <- 1L
  while (i <= length(args)) {
    a <- args[i]
    val <- if (i + 1L <= length(args)) args[i + 1L] else NA
    switch(a,
      "--out" = { out$out <- val; i <- i + 2L },
      stop("Unknown argument: ", a)
    )
  }
  out
}


load_csv_dir <- function(dir) {
  files <- list.files(dir, pattern = "\\.csv$", full.names = TRUE)
  if (length(files) == 0L) stop("No CSVs found in ", dir)
  do.call(rbind, lapply(files, utils::read.csv, stringsAsFactors = FALSE))
}


main <- function() {
  opts <- parse_cli_args()

  reps  <- load_csv_dir(RAW_DIR)
  truth <- load_csv_dir(TRUTH_DIR)

  truth_long <- bind_rows(
    truth %>% transmute(scenario_id, dgm, misspec, rho_L,
                        t, target = "S0", truth = S0),
    truth %>% transmute(scenario_id, dgm, misspec, rho_L,
                        t, target = "S1", truth = S1),
    truth %>% transmute(scenario_id, dgm, misspec, rho_L,
                        t, target = "RD", truth = RD)
  )

  joined <- reps %>%
    dplyr::left_join(
      truth_long %>% dplyr::select(scenario_id, t, target, truth),
      by = c("scenario_id", "t", "target")
    )

  summary <- joined %>%
    dplyr::group_by(scenario_id, dgm, misspec, rho_L, t, target) %>%
    dplyr::summarise(
      n_reps         = dplyr::n(),
      truth          = mean(truth, na.rm = TRUE),
      mean_est       = mean(est,   na.rm = TRUE),
      bias           = mean_est - truth,
      rel_bias_pct   = 100 * bias / abs(truth),
      empirical_sd   = stats::sd(est, na.rm = TRUE),
      mean_model_se  = mean(se, na.rm = TRUE),
      rel_se_err     = mean_model_se / empirical_sd - 1,
      coverage       = mean(ci_lo <= truth & truth <= ci_hi, na.rm = TRUE),
      power          = mean(!(ci_lo <= 0 & 0 <= ci_hi), na.rm = TRUE),
      mean_ci_width  = mean(ci_hi - ci_lo, na.rm = TRUE),
      .groups        = "drop"
    ) %>%
    dplyr::arrange(dgm, misspec, rho_L, target, t)

  utils::write.csv(summary, opts$out, row.names = FALSE)
  message(sprintf("[wrote] %s  (%d rows)", opts$out, nrow(summary)))

  # Console preview
  message("\n--- summary head ---")
  print(utils::head(summary, 20))
}

if (sys.nframe() == 0L) {
  main()
}
