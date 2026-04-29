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
# load_config.R checks for `yaml` only when read_dgm_yaml() is called,
# so sourcing it is harmless for PLR-only runs that never use --config.
source(file.path(CODE_DIR, "utils", "load_config.R"))
# AIPTW-Cox arm sourced lazily inside main() so a PLR-only run does
# not require riskRegression to be installed.


## ---- CLI argument parser ---------------------------------------------------

parse_cli_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  defaults <- list(
    R            = 200L,           # replications per scenario
    B            = 100L,           # bootstrap resamples per replicate
    base_seed    = 1000L,          # base seed; replicate r uses base_seed + r
    N            = 2500L,          # protocol Section 6.2
    t_eval       = c(1, 5, 10),    # protocol Section 6.3
    rescale_time = 0.25,           # discrete-time interval width; smaller = more time-period dummies in PLR (slower fit, finer grid)
    dgm          = NULL,           # NULL = all
    misspec      = NULL,           # NULL = all 5 main patterns
    rho          = NULL,           # NULL = all (0, 0.25, 0.75)
    include_heavy = FALSE,
    n_workers    = 1L,             # 1 = sequential
    overwrite    = FALSE,
    method       = "aiptw_plr",    # "aiptw_plr", "aiptw_cox", or "both"
    config       = NULL            # optional YAML config path
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
      "--rescale-time"  = { defaults$rescale_time <- as.numeric(val); i <- i + 2L },
      "--dgm"           = { defaults$dgm       <- strsplit(val, ",")[[1]];     i <- i + 2L },
      "--misspec"       = { defaults$misspec   <- strsplit(val, ",")[[1]];     i <- i + 2L },
      "--rho"           = { defaults$rho       <- as.numeric(strsplit(val, ",")[[1]]); i <- i + 2L },
      "--workers"       = { defaults$n_workers <- as.integer(val); i <- i + 2L },
      "--include-heavy" = { defaults$include_heavy <- TRUE; i <- i + 1L },
      "--overwrite"     = { defaults$overwrite <- TRUE; i <- i + 1L },
      "--method"        = { defaults$method <- val; i <- i + 2L },
      "--config"        = { defaults$config <- val; i <- i + 2L },
      { stop("Unknown CLI argument: ", a) }
    )
  }

  # Fail fast on obviously invalid --rescale-time. Per-DGM check (does
  # admin.cens divide evenly) happens at scenario start where admin.cens
  # is known.
  if (!is.finite(defaults$rescale_time) || defaults$rescale_time <= 0) {
    stop("--rescale-time must be a positive finite number; got ",
         defaults$rescale_time)
  }

  # Validate --method.
  if (!defaults$method %in% c("aiptw_plr", "aiptw_cox", "both")) {
    stop("--method must be one of 'aiptw_plr', 'aiptw_cox', 'both'; got ",
         defaults$method)
  }

  # --config can come in two flavours:
  #   (a) a SINGLE-DGM YAML (e.g. config/ph.yml). Then --dgm must be
  #       passed and must select exactly one DGM, otherwise the same
  #       PH parameter set would be applied to delayed/waning rows
  #       and produce mislabelled output.
  #   (b) per-DGM defaults at config/<dgm>.yml: not passed via
  #       --config; the loop picks them up automatically based on
  #       sc$dgm, so they apply correctly to a multi-DGM grid.
  if (!is.null(defaults$config)) {
    if (!file.exists(defaults$config)) {
      stop("--config file not found: ", defaults$config)
    }
    if (is.null(defaults$dgm) || length(defaults$dgm) != 1L) {
      stop("--config <yml> requires --dgm to be passed with exactly one ",
           "DGM (otherwise the YAML would be applied across mismatched ",
           "DGM labels). Pass --dgm ph (or delayed/waning), or omit ",
           "--config and let the per-DGM defaults at config/<dgm>.yml ",
           "apply automatically.")
    }
  }

  defaults
}


## Apply YAML-supplied DGM-level options on top of CLI args. Per-CLI
## overrides win: if the user passed --N or --t-eval explicitly we keep
## those, otherwise the YAML's values become the active values.
apply_yaml_to_opts <- function(opts, dgm_params, raw_argv) {
  cli_set <- function(flag) any(flag == raw_argv)
  if (!is.null(dgm_params$N) && !cli_set("--N")) {
    opts$N <- as.integer(dgm_params$N)
  }
  if (!is.null(dgm_params$t_eval) && !cli_set("--t-eval")) {
    opts$t_eval <- as.numeric(dgm_params$t_eval)
  }
  # N_truth has no CLI flag; YAML value wins if present.
  if (!is.null(dgm_params$N_truth)) {
    opts$N_truth <- as.integer(dgm_params$N_truth)
  }
  opts
}


## Methods to run for this invocation.
methods_to_run <- function(method_arg) {
  if (method_arg == "both") c("aiptw_plr", "aiptw_cox") else method_arg
}


## Lazy-load the AIPTW-Cox estimator only when needed (riskRegression is a
## heavy dependency we do not want to require for a PLR-only run).
ensure_cox_loaded <- function() {
  if (!exists("aiptw_cox_estimate", mode = "function")) {
    source(file.path(CODE_DIR, "utils", "aiptw_cox_estimator.R"))
  }
}


## ---- Parallel backend ------------------------------------------------------

set_parallel_plan <- function(n_workers) {
  if (n_workers <= 1L) {
    message("[parallel] sequential")
    return(invisible(1L))
  }
  has_pkgs <- requireNamespace("future", quietly = TRUE) &&
              requireNamespace("furrr",  quietly = TRUE)
  if (!has_pkgs) {
    warning("future / furrr not installed; falling back to sequential",
            call. = FALSE)
    return(invisible(1L))
  }
  # Respect container / cgroup CPU limits: parallelly::availableCores()
  # honours cgroups2.cpu.max so we do not over-subscribe a 1-CPU pod and
  # then fail in checkNumberOfLocalWorkers(). Cap the requested worker
  # count at the available cores with a clear warning.
  avail <- if (requireNamespace("parallelly", quietly = TRUE))
             as.integer(parallelly::availableCores()) else
             parallel::detectCores(logical = TRUE)
  if (is.na(avail) || avail < 1L) avail <- 1L
  if (n_workers > avail) {
    warning(sprintf(
      "--workers %d requested but only %d core(s) available; capping to %d",
      n_workers, avail, avail), call. = FALSE)
    n_workers <- avail
  }
  if (n_workers <= 1L) {
    message("[parallel] sequential (only ", avail, " core available)")
    return(invisible(1L))
  }
  future::plan(future::multisession, workers = n_workers)
  message("[parallel] future::multisession with ", n_workers, " workers")
  invisible(n_workers)
}

## Apply a function across replications either via furrr or a sequential
## lapply, returning a list. f(r) -> data.frame.
##
## The packages list is passed to furrr_options() so each multisession
## worker attaches them at startup. Without this, the rep_fn closure
## references unqualified symbols (mvrnorm, %>%, mutate, left_join,
## survSplit, ...) that are pulled in via source() in the parent session
## but not automatically reproduced on the worker, which would error or
## silently fall back to the wrong namespace.
##
## Sequential mode prints periodic progress messages so a long run on a
## 1-CPU pod does not appear frozen between scenario start and scenario
## complete.
run_replications <- function(R, f, n_workers = 1L, base_seed = 0L,
                              packages = c("survival", "dplyr", "MASS")) {
  if (n_workers > 1L && requireNamespace("furrr", quietly = TRUE)) {
    furrr::future_map(seq_len(R), f,
                      .options = furrr::furrr_options(
                        seed     = TRUE,
                        packages = packages
                      ))
  } else {
    t0       <- Sys.time()
    interval <- max(1L, R %/% 20L)         # ~20 progress messages per scenario
    out      <- vector("list", R)
    for (r in seq_len(R)) {
      out[[r]] <- f(r)
      if (r == 1L || r %% interval == 0L || r == R) {
        elapsed <- as.numeric(Sys.time() - t0, units = "secs")
        eta     <- elapsed * (R - r) / r
        message(sprintf(
          "    [rep %d/%d]  elapsed=%.1fs  eta=%.1fs", r, R, elapsed, eta
        ))
      }
    }
    out
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

run_one_scenario <- function(scenario_row, dgm_params, opts,
                              truth_cache = new.env(parent = emptyenv()),
                              rep_dir = NULL) {

  # Validate that admin.cens is an integer multiple of rescale_time. If
  # not, the discrete-time grid does not land on admin.cens and the
  # final evaluation point uses an earlier interval's Q while the IPCW
  # outcome is taken at admin.cens, biasing the t = admin.cens estimate.
  ratio <- dgm_params$admin.cens / opts$rescale_time
  if (abs(ratio - round(ratio)) > 1e-9) {
    stop(sprintf(
      "--rescale-time (%g) does not divide admin.cens (%g) evenly. Use a value such as 0.25, 0.5, 1, or 2.",
      opts$rescale_time, dgm_params$admin.cens
    ))
  }

  rho_L <- scenario_row$rho_L
  sigma_use <- make_sigma(rho_L = rho_L,
                           n_L = dgm_params$N.Lcovs.linear +
                                  dgm_params$N.Lcovs.sq)
  dgm_params$sigma <- sigma_use

  # Per-scenario rep-checkpoint directory. Caller supplies the full path
  # (which is method-subdir-aware); fall back to the legacy location for
  # in-process callers.
  if (is.null(rep_dir)) {
    rep_dir <- file.path(RAW_DIR, scenario_row$scenario_id)
  }
  dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)

  # Truth depends on (dgm, rho_L) but NOT on the misspecification
  # pattern, so cache by that key. Without caching, the 500k-row
  # truth simulation runs five times per (dgm, rho) pair (once per
  # misspec) which is wasted work.
  truth_key <- paste(scenario_row$dgm, scenario_row$rho_L, sep = "__")
  if (!exists(truth_key, envir = truth_cache, inherits = FALSE)) {
    base_truth <- compute_truth(
      t_eval     = opts$t_eval,
      dgm        = dgm_params,
      truth_seed = opts$base_seed * 31L + 1L
    )
    assign(truth_key, base_truth, envir = truth_cache)
  }
  truth <- get(truth_key, envir = truth_cache)
  # Annotate per-scenario keys (truth values themselves are identical
  # across misspec patterns within the same dgm/rho).
  truth$scenario_id <- scenario_row$scenario_id
  truth$dgm         <- scenario_row$dgm
  truth$misspec     <- scenario_row$misspec
  truth$rho_L       <- scenario_row$rho_L

  message(sprintf("[%s] truth: S0(t=10)=%.4f S1(t=10)=%.4f RD(t=10)=%.4f",
                  scenario_row$scenario_id,
                  truth$S0[truth$t == max(opts$t_eval)],
                  truth$S1[truth$t == max(opts$t_eval)],
                  truth$RD[truth$t == max(opts$t_eval)]))

  # Per-replication function (closure over opts, dgm_params, scenario).
  #
  # Runs each method requested by --method on the SAME simulated dataset
  # (so PLR and Cox arms share Monte-Carlo noise -- their differences
  # reflect only estimator choice, not data variation). Output rows are
  # tagged with a method column.
  active_methods <- methods_to_run(opts$method)

  rep_fn <- function(r) {
    rep_path <- file.path(rep_dir, sprintf("rep_%05d.csv", r))
    if (file.exists(rep_path)) {
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

    rep_rows <- vector("list", length(active_methods))

    for (mi in seq_along(active_methods)) {
      m <- active_methods[mi]

      if (m == "aiptw_plr") {
        est <- aiptw_estimate(
          surv.df, scenario_row$ps_spec, scenario_row$out_spec,
          t_eval = opts$t_eval, admin.cens = dgm_params$admin.cens,
          rescale_time = opts$rescale_time
        )
        boot <- aiptw_bootstrap(
          surv.df, scenario_row$ps_spec, scenario_row$out_spec,
          t_eval = opts$t_eval, admin.cens = dgm_params$admin.cens,
          rescale_time = opts$rescale_time, B = opts$B
        )
        # IF-based SE is not produced for the PLR arm; emit NA columns for
        # schema parity with the Cox arm.
        if_se_S0 <- NA_real_; if_se_S1 <- NA_real_; if_se_RD <- NA_real_
        if_lo_S0 <- NA_real_; if_hi_S0 <- NA_real_
        if_lo_S1 <- NA_real_; if_hi_S1 <- NA_real_
        if_lo_RD <- NA_real_; if_hi_RD <- NA_real_

      } else if (m == "aiptw_cox") {
        est <- aiptw_cox_estimate(
          surv.df, scenario_row$ps_spec, scenario_row$out_spec,
          t_eval = opts$t_eval, admin.cens = dgm_params$admin.cens
        )
        boot <- aiptw_cox_bootstrap(
          surv.df, scenario_row$ps_spec, scenario_row$out_spec,
          t_eval = opts$t_eval, admin.cens = dgm_params$admin.cens,
          B = opts$B
        )
      } else stop("unknown method: ", m)

      # Build the long-format rep_rows for this method. The per-row IF
      # columns come from the wide est data frame (Cox only); for PLR
      # they are NA-filled to keep one schema across both arms.
      method_long <- dplyr::bind_rows(
        data.frame(t = est$t, target = "S0", est = est$S0,
                   status = est$status, error_msg = est$error_msg,
                   if_se = if (m == "aiptw_cox") est$S0_if_se else NA_real_,
                   if_ci_lo = if (m == "aiptw_cox") est$S0_if_ci_lo else NA_real_,
                   if_ci_hi = if (m == "aiptw_cox") est$S0_if_ci_hi else NA_real_,
                   stringsAsFactors = FALSE),
        data.frame(t = est$t, target = "S1", est = est$S1,
                   status = est$status, error_msg = est$error_msg,
                   if_se = if (m == "aiptw_cox") est$S1_if_se else NA_real_,
                   if_ci_lo = if (m == "aiptw_cox") est$S1_if_ci_lo else NA_real_,
                   if_ci_hi = if (m == "aiptw_cox") est$S1_if_ci_hi else NA_real_,
                   stringsAsFactors = FALSE),
        data.frame(t = est$t, target = "RD", est = est$RD,
                   status = est$status, error_msg = est$error_msg,
                   if_se = if (m == "aiptw_cox") est$RD_if_se else NA_real_,
                   if_ci_lo = if (m == "aiptw_cox") est$RD_if_ci_lo else NA_real_,
                   if_ci_hi = if (m == "aiptw_cox") est$RD_if_ci_hi else NA_real_,
                   stringsAsFactors = FALSE)
      )
      method_long <- dplyr::left_join(method_long, boot,
                                       by = c("t", "target"))
      method_long$method      <- m
      method_long$rep         <- r
      method_long$scenario_id <- scenario_row$scenario_id
      method_long$dgm         <- scenario_row$dgm
      method_long$misspec     <- scenario_row$misspec
      method_long$rho_L       <- scenario_row$rho_L
      rep_rows[[mi]] <- method_long
    }

    rep_df <- do.call(rbind, rep_rows)

    # Atomic checkpoint: tmp -> rename. Checked because file.rename can
    # fail on Windows / networked storage.
    tmp_path <- paste0(rep_path, ".tmp")
    utils::write.csv(rep_df, tmp_path, row.names = FALSE)
    ok <- file.rename(tmp_path, rep_path)
    if (!isTRUE(ok)) {
      unlink(tmp_path)
      stop(sprintf("[checkpoint] file.rename failed: %s -> %s",
                   tmp_path, rep_path))
    }
    rep_df
  }

  rep_list <- run_replications(opts$R, rep_fn, opts$n_workers,
                                base_seed = opts$base_seed)

  rep_df <- do.call(rbind, rep_list)

  # Report failure rate up front. Count DISTINCT replicates that errored
  # (rep_df has 9 rows per replicate: 3 t x 3 target). Print the first
  # distinct error message inline so the user does not need to open the
  # raw CSV to diagnose a systematic failure.
  failed_reps <- unique(rep_df$rep[rep_df$status == "error"])
  n_fail <- length(failed_reps)
  if (n_fail > 0) {
    err_msgs <- unique(stats::na.omit(
      rep_df$error_msg[rep_df$status == "error"]
    ))
    first_err <- if (length(err_msgs) > 0L) err_msgs[1] else "(no message)"
    warning(sprintf(
      "[%s] %d / %d replicates failed. First error: %s",
      scenario_row$scenario_id, n_fail, opts$R, first_err
    ), call. = FALSE)
  }

  list(reps = rep_df, truth = truth, rep_dir = rep_dir)
}


## ---- Main ------------------------------------------------------------------

main <- function() {
  opts <- parse_cli_args()

  # Lazy-load the Cox arm dependencies if --method needs them.
  if (opts$method %in% c("aiptw_cox", "both")) ensure_cox_loaded()

  # Capture the actual worker count after parallelly capping.
  opts$n_workers <- set_parallel_plan(opts$n_workers)

  grid <- build_scenario_grid(
    dgms     = if (is.null(opts$dgm))     names(dgm_param_files) else opts$dgm,
    misspecs = if (is.null(opts$misspec)) names(misspec_patterns) else opts$misspec,
    rhos     = if (is.null(opts$rho))     rho_levels else opts$rho,
    include_heavy = opts$include_heavy
  )

  message(sprintf(
    "[grid] %d scenarios x R = %d reps x B = %d bootstraps  |  method = %s",
    nrow(grid), opts$R, opts$B, opts$method
  ))
  if (!is.null(opts$config)) {
    message("[config] loading DGM params from ", opts$config)
  }

  # DGM params are per-DGM (not per-rho); cache to avoid re-sourcing.
  # When --config is supplied, the YAML loader replaces the legacy
  # get_params_*.R lookup so all DGMs in --dgm draw from the same YAML
  # (or the per-DGM YAML files in config/).
  dgm_cache   <- new.env(parent = emptyenv())
  truth_cache <- new.env(parent = emptyenv())   # keyed by "<dgm>__<rho_L>"

  # Config that affects raw output content. Stored beside each scenario
  # CSV as a .cfg sidecar so we can detect and refuse to mix runs from
  # different settings (e.g. resuming a smoke run as a protocol run, or
  # changing --rescale-time / --N / --base-seed / --method mid-run).
  current_cfg <- list(
    N            = opts$N,
    R            = opts$R,
    B            = opts$B,
    rescale_time = opts$rescale_time,
    base_seed    = opts$base_seed,
    method       = opts$method,
    config       = if (is.null(opts$config)) "" else
                     normalizePath(opts$config, mustWork = FALSE)
  )

  # Method-specific subdirectory under results/raw/ so the PLR and Cox
  # arms can coexist without colliding on scenario_id. (Truth is method-
  # independent and stays under results/truth/.)
  method_subdir <- file.path(RAW_DIR, opts$method)
  dir.create(method_subdir, recursive = TRUE, showWarnings = FALSE)

  for (i in seq_len(nrow(grid))) {
    sc <- as.list(grid[i, ])

    raw_path   <- file.path(method_subdir, paste0(sc$scenario_id, ".csv"))
    truth_path <- file.path(TRUTH_DIR,     paste0(sc$scenario_id, ".csv"))
    rep_dir    <- file.path(method_subdir, sc$scenario_id)
    cfg_path   <- file.path(method_subdir, paste0(sc$scenario_id, ".cfg"))

    # Resolve DGM params FIRST (before the cfg check) so we can hash
    # their contents into the cfg. This catches YAML edits at the same
    # path -- previously the cfg only stored the path string, so an
    # in-place edit of config/ph.yml could be silently reused as if
    # compatible.
    if (!exists(sc$dgm, envir = dgm_cache, inherits = FALSE)) {
      # Three sources for the DGM parameter list, in priority:
      #   (1) --config <yml>: explicit path passed on the CLI
      #   (2) config/<dgm>.yml: per-DGM YAML in the config directory
      #   (3) get_params_<dgm>.R: legacy R parameter file (back-compat)
      yaml_default <- file.path(CODE_DIR, "config",
                                 paste0(sc$dgm, ".yml"))
      yaml_path <- if (!is.null(opts$config)) opts$config else
                    if (file.exists(yaml_default)) yaml_default else NULL

      if (!is.null(yaml_path)) {
        dp <- read_dgm_yaml(yaml_path)
        # YAML overrides for DGM-level options the CLI did not set
        # explicitly (--N, --t-eval, N_truth).
        opts <<- apply_yaml_to_opts(opts, dp,
                                     commandArgs(trailingOnly = TRUE))
        assign(sc$dgm, dp, envir = dgm_cache)
      } else {
        assign(sc$dgm, load_dgm_params(sc$dgm, code_dir = CODE_DIR),
               envir = dgm_cache)
      }
    }
    dgm_params <- get(sc$dgm, envir = dgm_cache)

    # Build the per-scenario config: run-level keys + a content-hash of
    # the DGM parameters so YAML edits are detected.
    scenario_cfg <- c(current_cfg,
                      list(dgm_hash = hash_dgm_params(dgm_params)))

    # Legacy checkpoint detection: outputs from a pre-.cfg run cannot
    # have their config validated. Refuse to resume without --overwrite.
    legacy_present <- (file.exists(raw_path) || dir.exists(rep_dir)) &&
                      !file.exists(cfg_path)
    if (legacy_present && !opts$overwrite) {
      stop(sprintf(
        "[%s] outputs exist from an older run with no .cfg sidecar; cannot validate config compatibility. Re-run with --overwrite or remove %s* first.",
        sc$scenario_id, raw_path
      ))
    }

    if (file.exists(cfg_path) && !opts$overwrite) {
      prev_cfg  <- dget(cfg_path)
      diff_keys <- names(scenario_cfg)[
        !mapply(identical, scenario_cfg,
                prev_cfg[names(scenario_cfg)])
      ]
      if (length(diff_keys) > 0L) {
        stop(sprintf(
          "[%s] existing scenario was produced with a different config (differs in: %s). Re-run with --overwrite or remove the stale outputs.",
          sc$scenario_id, paste(diff_keys, collapse = ", ")
        ))
      }
    }

    if (file.exists(raw_path) && !opts$overwrite) {
      message(sprintf("[skip] %s exists (use --overwrite to redo)", raw_path))
      next
    }

    # On --overwrite, clear any per-replicate checkpoint files too.
    if (opts$overwrite) {
      if (file.exists(raw_path))   unlink(raw_path)
      if (file.exists(truth_path)) unlink(truth_path)
      if (file.exists(cfg_path))   unlink(cfg_path)
      if (dir.exists(rep_dir))     unlink(rep_dir, recursive = TRUE)
    }

    # Write cfg before running so resume can validate.
    dput(scenario_cfg, file = cfg_path)

    message(sprintf("\n========== [%d/%d] %s ==========",
                    i, nrow(grid), sc$scenario_id))

    res <- run_one_scenario(sc, dgm_params, opts, truth_cache,
                             rep_dir = rep_dir)

    # Atomic consolidate. Order: truth first, then raw, because the
    # summariser uses raw existence as the "complete" signal and joins
    # to the truth file.
    tmp_truth <- paste0(truth_path, ".tmp")
    utils::write.csv(res$truth, tmp_truth, row.names = FALSE)
    ok_t <- file.rename(tmp_truth, truth_path)
    if (!isTRUE(ok_t)) {
      unlink(tmp_truth)
      stop(sprintf("[consolidate-truth] file.rename failed: %s -> %s",
                   tmp_truth, truth_path))
    }

    tmp_raw <- paste0(raw_path, ".tmp")
    utils::write.csv(res$reps, tmp_raw, row.names = FALSE)
    ok_r <- file.rename(tmp_raw, raw_path)
    if (!isTRUE(ok_r)) {
      unlink(tmp_raw)
      stop(sprintf("[consolidate-raw] file.rename failed: %s -> %s",
                   tmp_raw, raw_path))
    }

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
