## utils/scenarios.R
## Encodes the simulation scenario grid from Apr 28 methods document Section
## 6.4. A scenario is defined by:
##   - a DGM (proportional hazards / delayed effect / treatment waning)
##   - a misspecification pattern (which nuisance model is wrong, if any)
##   - an L-block correlation level (0 / 0.25 / 0.75)
##
## Each scenario also carries the PS and outcome covariate-spec strings to
## pass to aiptw_estimate(). These are aligned with the post-Apr-28 DGM:
##   PS truth   = linear in L1..L6 + W
##   Outcome    = linear L1..L6 + quadratic L1, L2 + O + W
## See file headers of sim_data.R and aiptw_discrete.R for the full mapping.

ps.specs <- c(
  correct  = "L1 + L2 + L3 + L4 + L5 + L6 + W",
  no_W     = "L1 + L2 + L3 + L4 + L5 + L6",
  wrong_ff = "L1sq + L2sq + L3 + L4 + L5 + L6 + W",
  heavy    = "L3 + L4 + L5"
)

out.specs <- c(
  correct  = "A + L1 + L1sq + L2 + L2sq + L3 + L4 + L5 + L6 + O + W",
  no_O     = "A + L1 + L1sq + L2 + L2sq + L3 + L4 + L5 + L6 + W",
  wrong_ff = "A + L1 + L2 + L3 + L4 + L5 + L6 + O + W",
  heavy    = "A + L3 + L4 + L5"
)

## Misspecification patterns: each isolates ONE nuisance model
misspec_patterns <- list(
  both_correct = list(ps = "correct",  out = "correct",
                       label = "Both correct"),
  miss_W_ps    = list(ps = "no_W",     out = "correct",
                       label = "PS: missing W"),
  miss_O_out   = list(ps = "correct",  out = "no_O",
                       label = "Outcome: missing O"),
  wrong_ff_out = list(ps = "correct",  out = "wrong_ff",
                       label = "Outcome: wrong form"),
  wrong_ff_ps  = list(ps = "wrong_ff", out = "correct",
                       label = "PS: wrong form")
)

## Correlation levels for the L block of Sigma
rho_levels <- c(none = 0, low = 0.25, high = 0.75)

## DGM scenarios: each maps to a get_params_*.R parameter file
dgm_param_files <- c(
  ph      = "get_params_ph.R",
  delayed = "get_params_delayed.R",
  waning  = "get_params_waning.R"
)


## Load a DGM parameter file into a self-contained list (no global pollution).
load_dgm_params <- function(dgm_name,
                             code_dir = ".") {
  stopifnot(dgm_name %in% names(dgm_param_files))
  env <- new.env(parent = baseenv())
  source(file.path(code_dir, dgm_param_files[[dgm_name]]),
         local = env, echo = FALSE)
  as.list(env)
}


## Build the full scenario grid as a tibble. Each row is one cell of the
## simulation study. Defaults match the protocol's main grid; pass
## include_heavy = TRUE to add the sensitivity heavy-misspec cells.
build_scenario_grid <- function(dgms       = c("ph", "delayed", "waning"),
                                  misspecs   = names(misspec_patterns),
                                  rhos       = rho_levels,
                                  include_heavy = FALSE) {

  if (include_heavy) {
    misspec_patterns_local <- c(
      misspec_patterns,
      list(
        heavy_ps  = list(ps = "heavy", out = "correct",
                          label = "PS: heavy misspec [sensitivity]"),
        heavy_out = list(ps = "correct", out = "heavy",
                          label = "Outcome: heavy misspec [sensitivity]")
      )
    )
  } else {
    misspec_patterns_local <- misspec_patterns
  }

  # Fail fast on unknown DGM / misspec names. silently dropping (e.g. via
  # intersect()) would let typos in --dgm / --misspec produce a smaller or
  # empty grid while the driver still reports "all scenarios complete".
  unknown_dgms <- setdiff(dgms, names(dgm_param_files))
  if (length(unknown_dgms) > 0L) {
    stop("Unknown dgm name(s): ", paste(unknown_dgms, collapse = ", "),
         ". Valid: ", paste(names(dgm_param_files), collapse = ", "))
  }
  unknown_misspecs <- setdiff(misspecs, names(misspec_patterns_local))
  if (length(unknown_misspecs) > 0L) {
    stop("Unknown misspec name(s): ",
         paste(unknown_misspecs, collapse = ", "),
         ". Valid: ", paste(names(misspec_patterns_local),
                             collapse = ", "))
  }
  if (length(misspecs) == 0L || length(dgms) == 0L || length(rhos) == 0L) {
    stop("Empty scenario grid (dgms / misspecs / rhos all required).")
  }

  rows <- expand.grid(
    dgm     = dgms,
    misspec = misspecs,
    rho_L   = rhos,
    stringsAsFactors = FALSE
  )

  rows$ps_spec  <- vapply(rows$misspec, function(m)
    ps.specs[[ misspec_patterns_local[[m]]$ps  ]], character(1))
  rows$out_spec <- vapply(rows$misspec, function(m)
    out.specs[[ misspec_patterns_local[[m]]$out ]], character(1))
  rows$label    <- vapply(rows$misspec, function(m)
    misspec_patterns_local[[m]]$label, character(1))

  rows$scenario_id <- sprintf("%s__%s__rho%g",
                               rows$dgm, rows$misspec, rows$rho_L)
  rows[, c("scenario_id", "dgm", "misspec", "rho_L",
           "ps_spec", "out_spec", "label")]
}
