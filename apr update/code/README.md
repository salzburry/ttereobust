# AIPTW simulation harness

Implementation of Section 6 of the Apr 28 methods document
(`apr update/Updated document Apr 28.pdf`) for the AIPTW arm. Two
estimator variants are produced side by side from the same simulated
data:

- **AIPTW-PLR** — discrete-time pooled-logistic outcome model with
  canonical logit link (`utils/aiptw_estimator.R`).
- **AIPTW-Cox** — Cox PH outcome model wrapping
  `riskRegression::ate()` (`utils/aiptw_cox_estimator.R`).

## File layout

```
apr update/code/
├── config/                  YAML DGM configs (option B, designed around
│   ├── default.yml          sim_surv_data() args)
│   ├── ph.yml
│   ├── delayed.yml
│   └── waning.yml
├── utils/
│   ├── sim_data.R           Section 6.2 data-generating mechanism
│   ├── cor_matrix.R         Sigma builder for the correlation sweep
│   ├── aiptw_estimator.R    PLR-AIPTW estimator + bootstrap
│   ├── aiptw_cox_estimator.R  Cox-AIPTW estimator + bootstrap
│   ├── compute_truth.R      Analytic Weibull g-formula truth
│   ├── scenarios.R          Scenario grid (DGM x misspec x rho_L)
│   └── load_config.R        YAML -> sim_surv_data params loader
├── get_params_ph.R          Legacy R parameter files (kept for
├── get_params_delayed.R     back-compat; YAML configs preferred)
├── get_params_waning.R
├── simulate_aiptw.R         Main replication driver
├── summarise_aiptw.R        Performance aggregator
├── run_aiptw_plr_analysis.R Wrapper script (PLR arm, protocol defaults)
├── run_aiptw_cox_analysis.R Wrapper script (Cox arm, protocol defaults)
├── run_pipeline.sh          End-to-end runner (simulate -> summarise -> render)
├── aiptw_discrete.R         Interactive demo: AIPTW-PLR across the 5
│                             main misspec scenarios on one dataset
├── aiptw.R                  Interactive demo: AIPTW-Cox single scenario
├── reports/
│   └── AIPTW_report.qmd     Quarto deliverable: HTML report consolidating
│                             tables + figures from results/summary.csv
└── results/                 (gitignored) raw replication output and summary
    ├── raw/
    ├── truth/
    └── summary.csv
```

## Dependencies

Install once on the pod:

```r
install.packages(c(
  # Core simulation + estimator
  "survival", "dplyr", "tidyr", "ggplot2", "MASS", "purrr",
  # YAML config support (--config flag)
  "yaml",
  # Cox AIPTW arm (riskRegression::ate())
  "riskRegression",
  # Optional: faster digest of DGM params for the cfg sidecar
  "digest",
  # Quarto report rendering
  "knitr", "rmarkdown"
))
```

Plus the `quarto` system binary (already present in most Domino R
images; otherwise install from <https://quarto.org/docs/get-started/>).

`yaml` and `riskRegression` are loaded only when the corresponding
features are exercised: a PLR-only run with no `--config` flag does
not require them.

## Run order

From `apr update/code/`:

```bash
# Single-fit timing test (one rep, no bootstrap, one scenario)
# Use this FIRST on any new pod to measure per-fit cost. Should
# finish in seconds. If not, the pod is CPU-throttled or
# --rescale-time needs to be coarsened.
Rscript simulate_aiptw.R --R 1 --B 0 \
  --dgm ph --misspec both_correct --rho 0 --overwrite

# Quick smoke test (a few scenarios, small R and B). Note this
# expands to 6 scenarios (2 misspec x 3 rho levels) x 50 reps x
# 51 fits per rep = 15,300 fits, so it is NOT a fast test on a
# single CPU. Add --rho 0 (single rho level) and/or
# --rescale-time 1 (10 time periods instead of 40) for a faster
# development run.
Rscript simulate_aiptw.R --R 50 --B 50 \
  --dgm ph --misspec both_correct,miss_W_ps

# Full protocol grid (3 DGMs x 5 misspec x 3 rho = 45 scenarios),
# AIPTW-PLR arm only
Rscript simulate_aiptw.R --R 1900 --B 200 --workers 8

# Both arms in one run (PLR + Cox on the same simulated datasets)
Rscript simulate_aiptw.R --method both --R 1900 --B 200 --workers 8

# YAML config (protocol-grade DGM from config/ph.yml)
Rscript simulate_aiptw.R --config config/ph.yml --method both \
                          --R 1900 --B 200 --workers 8

# Convention-aligned wrappers (one per arm)
Rscript run_aiptw_plr_analysis.R --workers 8
Rscript run_aiptw_cox_analysis.R --workers 8

# Aggregate into protocol metrics
Rscript summarise_aiptw.R

# End-to-end (simulate -> summarise -> render Quarto HTML report)
./run_pipeline.sh

# Render the Quarto report only (after simulate + summarise)
quarto render reports/AIPTW_report.qmd
```

### `--method` and `--config` flags

| Flag | Default | Effect |
|---|---|---|
| `--method aiptw_plr` | yes | discrete-time logit-PLR AIPTW (protocol estimator) |
| `--method aiptw_cox` | — | Cox PH outcome via `riskRegression::ate()` |
| `--method both` | — | run both arms on the same datasets; output rows tagged with `method` column |
| `--config <path>` | — | DGM parameters from a YAML config file. Overrides `get_params_*.R`. |

## Section -> implementation map

| Doc section | Implementation |
|---|---|
| 6.1 Aims | scope of the driver |
| 6.2 DGM (N = 2500, 8 covariates, exposure logit, Weibull hazard, indep censoring) | `utils/sim_data.R`, `utils/cor_matrix.R` |
| 6.3 Estimand (RD at t = 1, 5, 10) | `utils/aiptw_estimator.R` returns S0, S1, RD at `t_eval` |
| 6.4 Methods (AIPTW) and misspecification scenarios | `utils/aiptw_estimator.R` + `utils/scenarios.R` |
| 6.5 Performance (rel bias, coverage, rel SE error, power) | `summarise_aiptw.R` |
| Replications (1900) | `simulate_aiptw.R --R 1900` |

## CLI flags for `simulate_aiptw.R`

| Flag | Default | Purpose |
|---|---|---|
| `--R N` | 200 | replications per scenario |
| `--B N` | 100 | bootstrap resamples per replicate |
| `--N N` | 2500 | sample size per simulated dataset (Section 6.2) |
| `--rescale-time x` | 0.25 | discrete-time interval width in years; smaller = finer grid + more time-period dummies in the PLR (slower fit). Use 1.0 for ~4× faster runs on a constrained pod, 0.25 for protocol-faithful results |
| `--base-seed N` | 1000 | replication r uses `base_seed + r` |
| `--dgm a,b` | all | subset of `ph,delayed,waning` |
| `--misspec a,b` | all | subset of misspec patterns |
| `--rho a,b` | all | subset of `0,0.25,0.75` |
| `--workers N` | 1 | future::multisession workers (1 = sequential) |
| `--include-heavy` | off | add the sensitivity heavy-misspec cells |
| `--overwrite` | off | re-run scenarios whose CSV already exists |

## Output schema

`results/raw/<scenario_id>.csv` — one row per (replicate, t, target):

| column | description |
|---|---|
| scenario_id | `<dgm>__<misspec>__rho<rho_L>` |
| dgm, misspec, rho_L | scenario keys |
| method | `"aiptw_plr"` or `"aiptw_cox"` (which AIPTW variant produced this row) |
| rep | replicate index |
| t | evaluation time (1, 5, or 10) |
| target | `S0`, `S1`, or `RD` |
| est | point estimate (raw, unclipped) |
| status | `"ok"` or `"error"` |
| error_msg | error message when `status == "error"`, NA otherwise |
| n_boot_total | planned bootstrap reps (`B` from CLI) |
| n_boot_ok | successful bootstrap reps |
| n_boot_failed | `n_boot_total − n_boot_ok` |
| boot_errors | distinct error messages from any failed bootstrap fits, joined with ` | ` (NA when none) |
| se | bootstrap SE |
| ci_lo, ci_hi | 95% percentile bootstrap CI |
| if_se | influence-function SE (Cox arm only; NA on PLR rows) |
| if_ci_lo, if_ci_hi | 95% IF-based CI (Cox arm only; NA on PLR rows) |

`results/truth/<scenario_id>.csv` — one row per `t`, with `S0`, `S1`, and
`RD` as separate columns plus the scenario keys (`scenario_id`, `dgm`,
`misspec`, `rho_L`). The summariser pivots this into long form
(`t, target, truth`) before joining to the raw replicate output.

`results/summary.csv` — one row per (scenario, t, target) holding the
protocol Section 6.5 metrics plus audit columns:

| column | description |
|---|---|
| method | `"aiptw_plr"` or `"aiptw_cox"` (groups one row per variant) |
| n_reps_total | rows in the raw CSV for this group |
| n_reps_ok | replicates whose AIPTW estimator returned status = "ok" |
| n_reps_failed | `n_reps_total − n_reps_ok` |
| mean_n_boot_total | average planned bootstrap reps (= `B`) |
| mean_n_boot_ok | average successful bootstrap reps |
| mean_n_boot_failed | average failed bootstrap reps |
| n_ci_ok | replicates with finite bootstrap CI bounds (denominator for `coverage` and `power`) |
| n_if_ok | replicates with finite IF-based CI bounds (Cox arm only; denominator for `if_coverage`) |
| truth | analytic Weibull truth |
| mean_est | mean of `est` across replicates |
| bias, rel_bias_pct | absolute and relative bias |
| empirical_sd | sd of `est` across replicates |
| mean_model_se | mean of bootstrap SE |
| rel_se_err | `mean_model_se / empirical_sd − 1` |
| mean_if_se | mean of IF-based SE (Cox arm only; NA otherwise) |
| rel_if_se_err | `mean_if_se / empirical_sd − 1` (Cox arm only) |
| coverage | proportion of `n_ci_ok` bootstrap CIs containing the truth |
| if_coverage | proportion of `n_if_ok` IF-based CIs containing the truth (Cox arm only) |
| power | for `RD` rows only: proportion of CIs excluding 0; NA for `S0`, `S1` |
| mean_ci_width | average bootstrap CI width |
| boot_errors_seen | distinct bootstrap error messages observed across replicates, joined with ` || ` (NA when none) |

## Checkpointing and resume

Each replicate writes its row to `results/raw/<scenario_id>/rep_NNNNN.csv`
as soon as it completes. When all R replicates exist, the rep files are
consolidated into `results/raw/<scenario_id>.csv` (atomic via tmp+rename)
and the per-rep directory is removed.

If `simulate_aiptw.R` is interrupted partway through a scenario, restarting
it (without `--overwrite`) skips already-completed scenarios and resumes
the in-progress one from the last checkpointed replicate. Failed replicate
fits are recorded with `status = "error"` rather than dropped silently.
