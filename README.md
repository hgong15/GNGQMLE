# GNGQMLE replication materials

This repository contains data and code for the manuscript:

**Robust and Efficient Inference for GARCH Models with Heavy-Tailed and Asymmetric Innovations**

The package is organized to reproduce the simulation and empirical results reported in `main.tex` and `supp.tex`. Generated outputs are not stored in the repository; running the scripts will create `results/` and, where applicable, `figures/`.

## Directory structure

- `data/raw/`: raw public price/exchange-rate files used in the empirical applications.
- `data/processed/`: processed return series used in the manuscript tables.
- `R/`: scripts for preparing data, empirical applications, simulations, and figures/tables.
- `R/legacy/`: cleaned helper-function subsets, legacy source snapshots, and a legacy workspace required by the historical simulation/application code.

## Data

The main empirical application uses daily USD/TRY and USD/ARS exchange-rate data. The supplementary empirical application uses BTC, ETH, BNB, and TRX cryptocurrency price data. The raw files are included in `data/raw/`, and the processed return series used in the manuscript are included in `data/processed/`.

To regenerate the processed return files from the raw data, run from the repository root:

```r
Rscript R/create_processed_data.R
```

This writes:

- `data/processed/exchange_returns.csv`
- `data/processed/crypto_returns.csv`

The exchange-rate returns are computed as negative log price differences, matching the convention in the main empirical application. Cryptocurrency returns are computed as log price differences.

## Empirical application

Main exchange-rate application:

```r
Rscript R/final_exchange_application_results.R
Rscript R/compute_exchange_r_opt.R
```

Supplementary cryptocurrency application:

```r
Rscript R/final_crypto_application_diagnostics.R
```

These scripts use helper functions in `R/legacy/` and write CSV outputs under `results/application/`.

## Simulation results

### Estimator RMSE simulations

The legacy script used for the estimator-comparison simulations is included as:

- `R/legacy/Simulation_estimator_rmse.R`

This script is computationally intensive and was retained primarily for auditability of the historical simulation design.

### Main score-portmanteau simulation

To rerun the checkpointed simulation from scratch, use commands such as:

```r
Rscript R/run_main_portmanteau_r1_foreach_checkpoint.R --model=52 --replications=500
Rscript R/run_main_portmanteau_r1_foreach_checkpoint.R --model=53 --replications=500
```

These runs can be time-consuming. The scripts write checkpointed cell-level outputs and then combine them into `raw_results.csv` and `summary_stats.csv` under `results/simulation/portmanteau/`.

After the simulation summaries have been generated, publication-style PDFs can be created with commands such as:

```r
Rscript R/plot_main_portmanteau_publication.R --result_dir=results/simulation/portmanteau/main_portmanteau_r1_model52_rep500 --figure_prefix=test3_1
Rscript R/plot_main_portmanteau_publication.R --result_dir=results/simulation/portmanteau/main_portmanteau_r1_model53_rep500 --figure_prefix=test3_2
```

### Supplementary moment-based residual power diagnostics comparison

To rerun the supplementary comparison with moment-based residual power diagnostics, use:

```r
Rscript R/run_moment_power_comparison_foreach.R --n=500 --replications=1000
Rscript R/run_moment_power_comparison_foreach.R --n=1000 --replications=1000
Rscript R/combine_moment_power_comparison_tables.R --replications=1000
```

## R package requirements

The scripts use base R plus several CRAN packages, including:

- `foreach`
- `doParallel`
- `gamlss.dist`
- `numDeriv`
- `PearsonDS`
- `stabledist`

Some long-run simulation scripts use parallel workers. On Windows, they create PSOCK clusters through `doParallel`.

## Notes

The `R/legacy/` directory is included because parts of the simulation and application pipeline were originally developed using historical helper functions. The application scripts call the cleaned subset files `useful_fun2_subset.R` and `Eh1h2Er1r2_subset.R`; the original helper-source snapshots are kept as `.txt` files for auditability. The main replication scripts have been adjusted to use relative paths within this repository.
