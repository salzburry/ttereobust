## utils/load_config.R
##
## Reads a DGM YAML file and returns a list of parameters in the shape
## sim_surv_data() and compute_truth() already expect (legacy R parameter
## file names like coeff.A, coeff.L, alpha.L, gamma.tte, ...).
##
## The YAML schema is human-readable: hazard ratios + scaling factors are
## stored separately (mirroring the working group's config.yml shape) and
## the loader translates them into log-coefficients via
##     coeff.X = log(hazard_ratio.X) * scaling.X
##
## Only the `yaml` package is required at runtime. The loader tolerates
## a missing default.yml gracefully but will stop if the requested config
## omits required keys.

suppressPackageStartupMessages({
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required for --config support. ",
         "Install with install.packages('yaml').")
  }
})


## Deep-merge two named lists (override wins on key collision). Used to
## apply the per-DGM YAML on top of the shared defaults so any keys the
## per-DGM file omits inherit from default.yml.
deep_merge <- function(base, override) {
  if (is.null(override)) return(base)
  if (is.null(base))     return(override)
  for (k in names(override)) {
    if (is.list(base[[k]]) && is.list(override[[k]]) &&
        !is.null(names(base[[k]]))) {
      base[[k]] <- deep_merge(base[[k]], override[[k]])
    } else {
      base[[k]] <- override[[k]]
    }
  }
  base
}


## Read a YAML file and merge it on top of config/default.yml (if present).
read_dgm_yaml <- function(config_path,
                           default_path = file.path(
                             dirname(config_path), "default.yml"
                           )) {
  if (!file.exists(config_path)) {
    stop("Config file not found: ", config_path)
  }
  raw <- yaml::read_yaml(config_path)
  if (file.exists(default_path) &&
      normalizePath(default_path) != normalizePath(config_path)) {
    defaults <- yaml::read_yaml(default_path)
    raw <- deep_merge(defaults, raw)
  }
  validate_dgm_yaml(raw, config_path)
  yaml_to_dgm_params(raw)
}


## Validate that the parsed YAML has the keys we will read.
validate_dgm_yaml <- function(cfg, config_path) {
  required <- list(
    "dgm$N",                 cfg$dgm$N,
    "dgm$admin_cens",        cfg$dgm$admin_cens,
    "covariates$Lcovs_linear", cfg$covariates$Lcovs_linear,
    "covariates$Lcovs_sq",   cfg$covariates$Lcovs_sq,
    "covariates$mu",         cfg$covariates$mu,
    "exposure$alpha_L",      cfg$exposure$alpha_L,
    "exposure$alpha_W",      cfg$exposure$alpha_W,
    "outcome$hazard_ratios$A",   cfg$outcome$hazard_ratios$A,
    "outcome$hazard_ratios$L",   cfg$outcome$hazard_ratios$L,
    "outcome$hazard_ratios$Lsq", cfg$outcome$hazard_ratios$Lsq,
    "outcome$hazard_ratios$O",   cfg$outcome$hazard_ratios$O,
    "outcome$hazard_ratios$W",   cfg$outcome$hazard_ratios$W,
    "outcome$weibull$gamma",     cfg$outcome$weibull$gamma,
    "outcome$weibull$lambda",    cfg$outcome$weibull$lambda,
    "censoring$lambda_cens",     cfg$censoring$lambda_cens
  )
  keys  <- required[seq(1, length(required), 2)]
  vals  <- required[seq(2, length(required), 2)]
  missing <- vapply(vals, is.null, logical(1))
  if (any(missing)) {
    stop("Config '", config_path, "' is missing required keys: ",
         paste(unlist(keys[missing]), collapse = ", "))
  }
  invisible(TRUE)
}


## Translate the human-readable YAML into the legacy parameter list. A
## helper applies the optional per-block `scaling` to the log-hazard-ratio.
yaml_to_dgm_params <- function(cfg) {
  scl <- cfg$outcome$scaling %||% list()
  scale_log <- function(x, s) log(x) * (if (is.null(s)) 1 else s)

  k <- cfg$covariates$Lcovs_linear + cfg$covariates$Lcovs_sq + 2L
  mu <- as.numeric(cfg$covariates$mu)
  if (length(mu) != k) {
    stop("covariates$mu must have length ", k,
         " (Lcovs_linear + Lcovs_sq + 2); got ", length(mu))
  }

  list(
    # seed
    this.seed      = as.integer(cfg$seed %||% 2026),

    # DGM dimensions
    N.Lcovs.linear = cfg$covariates$Lcovs_linear,
    N.Lcovs.sq     = cfg$covariates$Lcovs_sq,
    mu             = mu,
    sigma          = diag(k),  # rho is applied later by make_sigma()
    admin.cens     = cfg$dgm$admin_cens,
    N              = cfg$dgm$N,
    N_truth        = cfg$dgm$N_truth %||% 500000L,

    # Exposure model
    alpha.L = as.numeric(cfg$exposure$alpha_L),
    alpha.W = as.numeric(cfg$exposure$alpha_W),

    # Outcome model log-HR coefficients (= log(HR) * scaling)
    coeff.A   = scale_log(cfg$outcome$hazard_ratios$A,   scl$A),
    coeff.L   = scale_log(cfg$outcome$hazard_ratios$L,   scl$L),
    coeff.Lsq = scale_log(cfg$outcome$hazard_ratios$Lsq, scl$Lsq),
    coeff.O   = scale_log(cfg$outcome$hazard_ratios$O,   scl$O),
    coeff.W   = scale_log(cfg$outcome$hazard_ratios$W,   scl$W),

    # Weibull baseline
    gamma.tte  = cfg$outcome$weibull$gamma,
    lambda.tte = cfg$outcome$weibull$lambda,

    # Censoring
    lambda.cens = cfg$censoring$lambda_cens,

    # Evaluation grid (overridable on the CLI)
    t_eval = as.numeric(cfg$evaluation$t_eval %||% c(1, 5, 10)),

    # Provenance
    config_path = normalizePath(attr(cfg, "config_path") %||%
                                 cfg$.config_path %||% "")
  )
}


## %||% — null-coalescing helper (rlang::%||% if available, else local).
`%||%` <- function(a, b) if (is.null(a)) b else a
