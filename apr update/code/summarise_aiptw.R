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
  # Read only the consolidated scenario files at the top of `dir`. Sub-
  # directories are per-rep checkpoints from in-progress simulate_aiptw.R
  # runs and should not be included in the aggregate.
  files <- list.files(dir, pattern = "\\.csv$", full.names = TRUE,
                       recursive = FALSE)
  if (length(files) == 0L) stop("No CSVs found in ", dir)
  do.call(rbind, lapply(files, utils::read.csv, stringsAsFactors = FALSE))
}


main <- function() {
  opts <- parse_cli_args()

  reps  <- load_csv_dir(RAW_DIR)
  truth <- load_csv_dir(TRUTH_DIR)

  # Some older outputs may not have status / n_boot_ok columns; default them.
  if (is.null(reps$status))    reps$status    <- "ok"
  if (is.null(reps$n_boot_ok)) reps$n_boot_ok <- NA_integer_

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
      n_reps_total   = dplyr::n(),
      n_reps_ok      = sum(status == "ok", na.rm = TRUE),
      n_reps_failed  = n_reps_total - n_reps_ok,
      mean_n_boot_ok = mean(n_boot_ok, na.rm = TRUE),
      truth          = mean(truth, na.rm = TRUE),
      mean_est       = mean(est,   na.rm = TRUE),
      bias           = mean_est - truth,
      rel_bias_pct   = 100 * bias / abs(truth),
      empirical_sd   = stats::sd(est, na.rm = TRUE),
      mean_model_se  = mean(se, na.rm = TRUE),
      rel_se_err     = mean_model_se / empirical_sd - 1,
      coverage       = mean(ci_lo <= truth & truth <= ci_hi, na.rm = TRUE),
      # Power: P(CI excludes 0). Meaningful only for the RD target where
      # H0: RD = 0 makes sense. S0 and S1 are survival probabilities far
      # from 0 by construction, so 'power' for those targets is reported
      # as NA to avoid misinterpretation.
      power          = ifelse(
        dplyr::first(target) == "RD",
        mean(!(ci_lo <= 0 & 0 <= ci_hi), na.rm = TRUE),
        NA_real_
      ),
      mean_ci_width  = mean(ci_hi - ci_lo, na.rm = TRUE),
      .groups        = "drop"
    ) %>%
    dplyr::arrange(dgm, misspec, rho_L, target, t)

  utils::write.csv(summary, opts$out, row.names = FALSE)
  message(sprintf("[wrote] %s  (%d rows)", opts$out, nrow(summary)))

  n_failed_total <- sum(summary$n_reps_failed[summary$target == "S0"])
  if (n_failed_total > 0) {
    warning(sprintf("%d total replicate failures across scenarios; ",
                    n_failed_total),
            "inspect status / error_msg in results/raw/<scenario>.csv",
            call. = FALSE)
  }

  # Console preview
  message("\n--- summary head ---")
  print(utils::head(summary, 20))
}

if (sys.nframe() == 0L) {
  main()
}
