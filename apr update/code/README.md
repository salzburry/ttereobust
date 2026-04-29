# AIPTW simulation harness

Implementation of Section 6 of the Apr 28 methods document
(`apr update/Updated document Apr 28.pdf`) for the AIPTW arm.

## File layout

```
apr update/code/
├── utils/
│   ├── sim_data.R           Section 6.2 data-generating mechanism
│   ├── cor_matrix.R         Sigma builder for the correlation sweep
│   ├── aiptw_estimator.R    Pure AIPTW estimator + bootstrap
│   ├── compute_truth.R      Analytic Weibull g-formula truth
│   └── scenarios.R          Scenario grid (DGM x misspec x rho_L)
├── get_params_ph.R          PH DGM parameters
├── get_params_delayed.R     Delayed-effect DGM parameters
├── get_params_waning.R      Waning-effect DGM parameters
├── simulate_aiptw.R         Main replication driver
├── summarise_aiptw.R        Performance aggregator
├── aiptw_discrete.R         Interactive demo: AIPTW-PLR across the 5
│                             main misspec scenarios on one dataset
│                             (uses utils/aiptw_estimator.R + scenarios.R).
├── aiptw.R                  Interactive demo: AIPTW-Cox (continuous-time
│                             outcome variant; not in the harness).
└── results/                 (gitignored) raw replication output and summary
    ├── raw/
    ├── truth/
    └── summary.csv
```

## Run order

From `apr update/code/`:

```bash
# Smoke test (a few scenarios, small R and B)
Rscript simulate_aiptw.R --R 50 --B 50 --dgm ph --misspec both_correct,miss_W_ps

# Full protocol grid (DGMs x misspec x rho)
Rscript simulate_aiptw.R --R 1900 --B 200 --workers 8

# Aggregate into protocol metrics
Rscript summarise_aiptw.R
```

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

`results/truth/<scenario_id>.csv` — one row per `t`, with `S0`, `S1`, and
`RD` as separate columns plus the scenario keys (`scenario_id`, `dgm`,
`misspec`, `rho_L`). The summariser pivots this into long form
(`t, target, truth`) before joining to the raw replicate output.

`results/summary.csv` — one row per (scenario, t, target) holding the
protocol Section 6.5 metrics plus audit columns:

| column | description |
|---|---|
| n_reps_total | rows in the raw CSV for this group |
| n_reps_ok | replicates whose AIPTW estimator returned status = "ok" |
| n_reps_failed | `n_reps_total − n_reps_ok` |
| mean_n_boot_total | average planned bootstrap reps (= `B`) |
| mean_n_boot_ok | average successful bootstrap reps |
| mean_n_boot_failed | average failed bootstrap reps |
| n_ci_ok | replicates with finite CI bounds (denominator for `coverage` and `power`) |
| truth | analytic Weibull truth |
| mean_est | mean of `est` across replicates |
| bias, rel_bias_pct | absolute and relative bias |
| empirical_sd | sd of `est` across replicates |
| mean_model_se | mean of bootstrap SE |
| rel_se_err | `mean_model_se / empirical_sd − 1` |
| coverage | proportion of `n_ci_ok` CIs containing the truth |
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
