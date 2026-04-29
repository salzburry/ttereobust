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
├── aiptw_discrete.R         Single-scenario interactive demo
├── aiptw.R                  Continuous-time / Cox-outcome AIPTW variant
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
| se | bootstrap SE |
| ci_lo, ci_hi | 95% percentile bootstrap CI |

`results/truth/<scenario_id>.csv` — one row per (t, target) holding the
analytic Weibull truth.

`results/summary.csv` — one row per (scenario, t, target) holding the
protocol Section 6.5 metrics.
