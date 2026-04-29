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

  # Some older outputs may not have the audit columns; default them so
  # the aggregator still produces a usable summary.
  if (is.null(reps$status))        reps$status        <- "ok"
  if (is.null(reps$n_boot_ok))     reps$n_boot_ok     <- NA_integer_
  if (is.null(reps$n_boot_total))  reps$n_boot_total  <- NA_integer_
  if (is.null(reps$n_boot_failed)) reps$n_boot_failed <- NA_integer_
  if (is.null(reps$boot_errors))   reps$boot_errors   <- NA_character_

  truth_long <- bind_rows(
    truth %>% transmute(scenario_id, dgm, misspec, rho_L,
                        t, target = "S0", truth = S0),
    truth %>% transmute(scenario_id, dgm, misspec, rho_L,
                        t, target = "S1", truth = S1),
    truth %>% transmute(scenario_id, dgm, misspec, rho_L,
                        t, target = "RD", truth = RD)
  )

  # Validate raw <-> truth alignment BEFORE the join so a stale or missing
  # truth file fails loudly instead of silently producing NA bias / NaN
  # coverage downstream.
  #
  # Two layers of validation:
  #   (1) Scenario-level: every raw scenario_id must appear in truth.
  #   (2) Row-level: every (scenario_id, t, target) triple in the raw
  #       output must match a row in truth_long. A truth CSV that has
  #       the right scenario_id but is missing t = 5 or t = 10 (e.g. an
  #       older partial output) would otherwise produce NA truth and
  #       NaN performance metrics after the join.
  raw_scenarios   <- unique(reps$scenario_id)
  truth_scenarios <- unique(truth$scenario_id)
  missing_truth   <- setdiff(raw_scenarios, truth_scenarios)
  if (length(missing_truth) > 0L) {
    stop("Raw scenarios with no matching truth CSV (would produce NA ",
         "performance metrics): ",
         paste(missing_truth, collapse = ", "),
         ". Re-run simulate_aiptw.R for these scenarios or remove the ",
         "stale raw outputs.")
  }
  orphan_truth <- setdiff(truth_scenarios, raw_scenarios)
  if (length(orphan_truth) > 0L) {
    warning("Truth CSVs with no matching raw scenario (will be ignored): ",
            paste(orphan_truth, collapse = ", "), call. = FALSE)
  }

  reps_keys  <- dplyr::distinct(
    reps[, c("scenario_id", "t", "target"), drop = FALSE]
  )
  truth_keys <- dplyr::distinct(
    truth_long[, c("scenario_id", "t", "target"), drop = FALSE]
  )
  missing_rows <- dplyr::anti_join(
    reps_keys, truth_keys,
    by = c("scenario_id", "t", "target")
  )
  if (nrow(missing_rows) > 0L) {
    msg <- utils::capture.output(print(utils::head(missing_rows, 10)))
    stop("Truth rows missing for some (scenario_id, t, target) triples ",
         "(stale or partial truth CSV). First missing rows:\n",
         paste(msg, collapse = "\n"))
  }

  joined <- reps %>%
    dplyr::left_join(
      truth_long %>% dplyr::select(scenario_id, t, target, truth),
      by = c("scenario_id", "t", "target")
    )

  summary <- joined %>%
    dplyr::group_by(scenario_id, dgm, misspec, rho_L, t, target) %>%
    dplyr::summarise(
      # Replicate-level accounting (denominator audit).
      n_reps_total      = dplyr::n(),
      n_reps_ok         = sum(status == "ok", na.rm = TRUE),
      n_reps_failed     = n_reps_total - n_reps_ok,
      # Bootstrap-level accounting (denominator for SE/CI).
      mean_n_boot_total = mean(n_boot_total, na.rm = TRUE),
      mean_n_boot_ok    = mean(n_boot_ok,    na.rm = TRUE),
      mean_n_boot_failed = mean(n_boot_failed, na.rm = TRUE),
      # CI denominator: replicates whose bootstrap actually produced a
      # finite CI. coverage and power below average over these rows; the
      # denominator is stored explicitly so it is auditable.
      n_ci_ok           = sum(!is.na(ci_lo) & !is.na(ci_hi)),
      truth             = mean(truth, na.rm = TRUE),
      mean_est          = mean(est,   na.rm = TRUE),
      bias              = mean_est - truth,
      rel_bias_pct      = 100 * bias / abs(truth),
      empirical_sd      = stats::sd(est, na.rm = TRUE),
      mean_model_se     = mean(se, na.rm = TRUE),
      rel_se_err        = mean_model_se / empirical_sd - 1,
      coverage          = mean(ci_lo <= truth & truth <= ci_hi,
                                na.rm = TRUE),
      # Power: P(CI excludes 0). Meaningful only for the RD target where
      # H0: RD = 0 makes sense. S0 and S1 are survival probabilities far
      # from 0 by construction, so 'power' for those targets is NA.
      power             = ifelse(
        dplyr::first(target) == "RD",
        mean(!(ci_lo <= 0 & 0 <= ci_hi), na.rm = TRUE),
        NA_real_
      ),
      mean_ci_width     = mean(ci_hi - ci_lo, na.rm = TRUE),
      # Distinct bootstrap error messages observed across replicates,
      # joined with " || " across reps. Lets a high n_reps_failed or
      # mean_n_boot_failed be diagnosed without re-running.
      boot_errors_seen  = {
        msgs <- unique(stats::na.omit(boot_errors))
        if (length(msgs) == 0L) NA_character_ else
          paste(msgs, collapse = " || ")
      },
      .groups           = "drop"
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
