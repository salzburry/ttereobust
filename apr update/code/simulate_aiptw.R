## simulate_aiptw.R
## Replication driver for the AIPTW arm of the Apr 28 simulation study.
##
## Implements the protocol-level loop described in Section 6 of the methods
## document:
##   - Section 6.2 DGM:                 utils/sim_data.R + utils/cor_matrix.R
##   - Section 6.3 Estimand (RD at 1, 5, 10):  utils/aiptw_estimator.R
##   - Section 6.4 Methods (AIPTW):     utils/aiptw_estimator.R + scenarios.R
##   - Section 6.5 Performance:         summarise_aiptw.R (separate script)
##
## For each scenario in the protocol grid (DGM x misspecification x rho_L):
##   1. Compute the analytic Weibull truth at t = 1, 5, 10.
##   2. For r = 1..R replications:
##         - Simulate one dataset of N = 2500.
##         - Run AIPTW for the scenario's PS and outcome specs.
##         - Run B bootstrap resamples to obtain SE and percentile CIs.
##         - Append (rep, t, target, est, se, ci_lo, ci_hi) to the scenario
##           result table.
##   3. Write the scenario result CSV under results/raw/.
##   4. Write the truth CSV under results/truth/.
##
## Defaults: R = 200, B = 100 (development-fast). Bump R to 1900 and B to 200
## for the protocol-grade run.
##
## Parallel execution: optional via {future}/{furrr}; falls back to sequential
## if those packages are not installed. Set N_WORKERS via Sys.setenv() or by
## editing the call to set_parallel_plan() below.
##
## CLI usage (from apr update/code/):
##   Rscript simulate_aiptw.R                        # all scenarios, defaults
##   Rscript simulate_aiptw.R --R 50 --B 50          # quick smoke test
##   Rscript simulate_aiptw.R --dgm ph --misspec both_correct,miss_W_ps
## (See parse_cli_args() below for full list of flags.)

suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
  library(MASS)
})

## ---- Configuration ---------------------------------------------------------

CODE_DIR    <- getwd()                   # run from apr update/code/
RESULTS_DIR <- file.path(CODE_DIR, "results")
RAW_DIR     <- file.path(RESULTS_DIR, "raw")
TRUTH_DIR   <- file.path(RESULTS_DIR, "truth")

dir.create(RAW_DIR,   recursive = TRUE, showWarnings = FALSE)
dir.create(TRUTH_DIR, recursive = TRUE, showWarnings = FALSE)

source(file.path(CODE_DIR, "utils", "sim_data.R"))
source(file.path(CODE_DIR, "utils", "cor_matrix.R"))
source(file.path(CODE_DIR, "utils", "aiptw_estimator.R"))
source(file.path(CODE_DIR, "utils", "compute_truth.R"))
source(file.path(CODE_DIR, "utils", "scenarios.R"))


## ---- CLI argument parser ---------------------------------------------------

parse_cli_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  defaults <- list(
    R         = 200L,            # replications per scenario
    B         = 100L,            # bootstrap resamples per replicate
    base_seed = 1000L,           # base seed; replicate r uses base_seed + r
    N         = 2500L,           # protocol Section 6.2
    t_eval    = c(1, 5, 10),     # protocol Section 6.3
    dgm       = NULL,            # NULL = all
    misspec   = NULL,            # NULL = all 5 main patterns
    rho       = NULL,            # NULL = all (0, 0.25, 0.75)
    include_heavy = FALSE,
    n_workers = 1L,              # 1 = sequential
    overwrite = FALSE
  )

  i <- 1L
  while (i <= length(args)) {
    a <- args[i]
    val <- if (i + 1L <= length(args)) args[i + 1L] else NA
    switch(a,
      "--R"             = { defaults$R         <- as.integer(val); i <- i + 2L },
      "--B"             = { defaults$B         <- as.integer(val); i <- i + 2L },
      "--base-seed"     = { defaults$base_seed <- as.integer(val); i <- i + 2L },
      "--N"             = { defaults$N         <- as.integer(val); i <- i + 2L },
      "--dgm"           = { defaults$dgm       <- strsplit(val, ",")[[1]];     i <- i + 2L },
      "--misspec"       = { defaults$misspec   <- strsplit(val, ",")[[1]];     i <- i + 2L },
      "--rho"           = { defaults$rho       <- as.numeric(strsplit(val, ",")[[1]]); i <- i + 2L },
      "--workers"       = { defaults$n_workers <- as.integer(val); i <- i + 2L },
      "--include-heavy" = { defaults$include_heavy <- TRUE; i <- i + 1L },
      "--overwrite"     = { defaults$overwrite <- TRUE; i <- i + 1L },
      { stop("Unknown CLI argument: ", a) }
    )
  }
  defaults
}


## ---- Parallel backend ------------------------------------------------------

set_parallel_plan <- function(n_workers) {
  if (n_workers <= 1L) {
    message("[parallel] sequential")
    return(invisible(NULL))
  }
  has_pkgs <- requireNamespace("future", quietly = TRUE) &&
              requireNamespace("furrr",  quietly = TRUE)
  if (!has_pkgs) {
    warning("future / furrr not installed; falling back to sequential",
            call. = FALSE)
    return(invisible(NULL))
  }
  future::plan(future::multisession, workers = n_workers)
  message("[parallel] future::multisession with ", n_workers, " workers")
}

## Apply a function across replications either via furrr or a sequential
## lapply, returning a list. f(r) -> data.frame.
run_replications <- function(R, f, n_workers = 1L, base_seed = 0L) {
  if (n_workers > 1L && requireNamespace("furrr", quietly = TRUE)) {
    furrr::future_map(seq_len(R), f,
                      .options = furrr::furrr_options(seed = TRUE))
  } else {
    lapply(seq_len(R), f)
  }
}


## ---- Single-scenario runner ------------------------------------------------
##
## Crash-safe checkpoint protocol:
##   - Each replicate writes its result to results/raw/<scenario_id>/rep_NNNNN.csv
##     immediately on completion. If any rep file already exists when the
##     scenario starts, it is loaded as-is and the rep is not re-run; this
##     enables resume-after-interrupt without --overwrite.
##   - When all R replicates have files on disk, the rep files are
##     consolidated into a single results/raw/<scenario_id>.csv and the
##     per-rep directory is removed. The summariser reads only the
##     consolidated CSVs, so an interrupted scenario (per-rep dir present,
##     consolidated CSV absent) is correctly skipped by the aggregator until
##     the run completes.

run_one_scenario <- function(scenario_row, dgm_params, opts) {

  rho_L <- scenario_row$rho_L
  sigma_use <- make_sigma(rho_L = rho_L,
                           n_L = dgm_params$N.Lcovs.linear +
                                  dgm_params$N.Lcovs.sq)
  dgm_params$sigma <- sigma_use

  # Per-scenario rep-checkpoint directory
  rep_dir <- file.path(RAW_DIR, scenario_row$scenario_id)
  dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)

  # Truth (does not depend on rep seed)
  truth <- compute_truth(
    t_eval     = opts$t_eval,
    dgm        = dgm_params,
    truth_seed = opts$base_seed * 31L + 1L
  )
  truth$scenario_id <- scenario_row$scenario_id
  truth$dgm         <- scenario_row$dgm
  truth$misspec     <- scenario_row$misspec
  truth$rho_L       <- scenario_row$rho_L

  message(sprintf("[%s] truth: S0(t=10)=%.4f S1(t=10)=%.4f RD(t=10)=%.4f",
                  scenario_row$scenario_id,
                  truth$S0[truth$t == max(opts$t_eval)],
                  truth$S1[truth$t == max(opts$t_eval)],
                  truth$RD[truth$t == max(opts$t_eval)]))

  # Per-replication function (closure over opts and scenario)
  rep_fn <- function(r) {
    rep_path <- file.path(rep_dir, sprintf("rep_%05d.csv", r))
    if (file.exists(rep_path)) {
      # Resume: rep already completed in a previous invocation.
      return(utils::read.csv(rep_path, stringsAsFactors = FALSE))
    }

    sim <- sim_surv_data(
      seed         = opts$base_seed + r,
      N            = opts$N,
      Lcovs.linear = dgm_params$N.Lcovs.linear,
      Lcovs.sq     = dgm_params$N.Lcovs.sq,
      mu           = dgm_params$mu,
      sigma        = dgm_params$sigma,
      alpha.L      = dgm_params$alpha.L,
      alpha.W      = dgm_params$alpha.W,
      coeff.A      = dgm_params$coeff.A,
      coeff.L      = dgm_params$coeff.L,
      coeff.Lsq    = dgm_params$coeff.Lsq,
      coeff.O      = dgm_params$coeff.O,
      coeff.W      = dgm_params$coeff.W,
      gamma.tte    = dgm_params$gamma.tte,
      lambda.tte   = dgm_params$lambda.tte,
      lambda.cens  = dgm_params$lambda.cens,
      admin.cens   = dgm_params$admin.cens,
      gen.truth    = NA
    )
    surv.df <- sim$data

    est <- aiptw_estimate(
      surv.df, scenario_row$ps_spec, scenario_row$out_spec,
      t_eval     = opts$t_eval,
      admin.cens = dgm_params$admin.cens
    )

    boot <- aiptw_bootstrap(
      surv.df, scenario_row$ps_spec, scenario_row$out_spec,
      t_eval     = opts$t_eval,
      admin.cens = dgm_params$admin.cens,
      B          = opts$B
    )

    # Long-format with status carried through.
    est_long <- dplyr::bind_rows(
      data.frame(t = est$t, target = "S0", est = est$S0,
                 status = est$status, error_msg = est$error_msg,
                 stringsAsFactors = FALSE),
      data.frame(t = est$t, target = "S1", est = est$S1,
                 status = est$status, error_msg = est$error_msg,
                 stringsAsFactors = FALSE),
      data.frame(t = est$t, target = "RD", est = est$RD,
                 status = est$status, error_msg = est$error_msg,
                 stringsAsFactors = FALSE)
    )

    rep_df <- dplyr::left_join(est_long, boot, by = c("t", "target"))
    rep_df$rep         <- r
    rep_df$scenario_id <- scenario_row$scenario_id
    rep_df$dgm         <- scenario_row$dgm
    rep_df$misspec     <- scenario_row$misspec
    rep_df$rho_L       <- scenario_row$rho_L

    # Atomic checkpoint: write to a tmp file, then rename.
    tmp_path <- paste0(rep_path, ".tmp")
    utils::write.csv(rep_df, tmp_path, row.names = FALSE)
    file.rename(tmp_path, rep_path)
    rep_df
  }

  rep_list <- run_replications(opts$R, rep_fn, opts$n_workers,
                                base_seed = opts$base_seed)

  rep_df <- do.call(rbind, rep_list)

  # Report failure rate up front so the user notices early.
  n_fail <- sum(rep_df$status == "error" & rep_df$target == "S0")
  if (n_fail > 0) {
    warning(sprintf("[%s] %d / %d replicates failed; see error_msg column",
                    scenario_row$scenario_id, n_fail, opts$R),
            call. = FALSE)
  }

  list(reps = rep_df, truth = truth, rep_dir = rep_dir)
}


## ---- Main ------------------------------------------------------------------

main <- function() {
  opts <- parse_cli_args()
  set_parallel_plan(opts$n_workers)

  grid <- build_scenario_grid(
    dgms     = if (is.null(opts$dgm))     names(dgm_param_files) else opts$dgm,
    misspecs = if (is.null(opts$misspec)) names(misspec_patterns) else opts$misspec,
    rhos     = if (is.null(opts$rho))     rho_levels else opts$rho,
    include_heavy = opts$include_heavy
  )

  message(sprintf("[grid] %d scenarios x R = %d reps x B = %d bootstraps",
                  nrow(grid), opts$R, opts$B))

  # DGM params are per-DGM (not per-rho); cache to avoid re-sourcing
  dgm_cache <- new.env(parent = emptyenv())

  for (i in seq_len(nrow(grid))) {
    sc <- as.list(grid[i, ])

    raw_path   <- file.path(RAW_DIR,   paste0(sc$scenario_id, ".csv"))
    truth_path <- file.path(TRUTH_DIR, paste0(sc$scenario_id, ".csv"))
    if (file.exists(raw_path) && !opts$overwrite) {
      message(sprintf("[skip] %s exists (use --overwrite to redo)", raw_path))
      next
    }

    if (!exists(sc$dgm, envir = dgm_cache, inherits = FALSE)) {
      assign(sc$dgm, load_dgm_params(sc$dgm, code_dir = CODE_DIR),
             envir = dgm_cache)
    }
    dgm_params <- get(sc$dgm, envir = dgm_cache)

    message(sprintf("\n========== [%d/%d] %s ==========",
                    i, nrow(grid), sc$scenario_id))

    res <- run_one_scenario(sc, dgm_params, opts)

    # Atomic consolidate: write tmp -> rename, then drop the rep_dir so
    # the summariser sees only completed scenarios.
    tmp_raw <- paste0(raw_path, ".tmp")
    utils::write.csv(res$reps,  tmp_raw,    row.names = FALSE)
    file.rename(tmp_raw, raw_path)

    utils::write.csv(res$truth, truth_path, row.names = FALSE)
    if (dir.exists(res$rep_dir)) {
      unlink(res$rep_dir, recursive = TRUE)
    }
    message(sprintf("[wrote] %s  (%d rows)", raw_path, nrow(res$reps)))
  }

  message("\n[done] all scenarios complete")
}

if (sys.nframe() == 0L) {
  main()
}
