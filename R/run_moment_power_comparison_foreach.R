## Compact comparison with moment-based residual power diagnostics.
##
## Design:
##   - DGP: model 5.2 in the main text
##   - Estimator: the same old-workspace GNGQMLE residuals, r = 1
##   - Tests: proposed score test vs. a moment-based residual power diagnostic with power s = 1
##   - n = 500, 200 replications
##   - innovations: N(0,1) and t3
##   - baselines: stationary, boundary, explosive by default
##   - z0: 0, 0.3, 0.6 by default
##
## The script uses foreach + doParallel.  On Windows this creates a PSOCK
## cluster; each worker loads the historical .RData workspace directly.

suppressPackageStartupMessages({
  library(foreach)
  library(doParallel)
})

old_workspace <- "R/legacy/legacy_workspace.RData"
args <- commandArgs(trailingOnly = TRUE)

get_arg_value <- function(name, default = NA_character_) {
  pattern <- paste0("^--", name, "=")
  matched <- grep(pattern, args, value = TRUE)
  if (length(matched) == 0L) {
    return(default)
  }
  sub(pattern, "", matched[1])
}

replications <- as.integer(get_arg_value("replications", "200"))
n_samples <- as.integer(get_arg_value("n", "500"))
parse_csv_arg <- function(name, default) {
  value <- get_arg_value(name, default)
  trimws(strsplit(value, ",", fixed = TRUE)[[1]])
}

out_dir <- get_arg_value(
  "outdir",
  sprintf("results/simulation/moment_power_comparison/moment_power_comparison_foreach_r1_n%d_rep%d", n_samples, replications)
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

m_lag <- 6L
r_power <- 1
s_power <- 1
ngqmle_params <- c(1, 7)
z_grid <- as.numeric(parse_csv_arg("z_grid", "0,0.3,0.6"))
regime_grid <- parse_csv_arg("regimes", "stationary,boundary,explosive")
seed_base <- 20260419L

innovation_specs <- list(
  norm = list(
    data_params = c(0, 2),
    alpha = c(stationary = 0.06, boundary = 0.1096508, explosive = 0.16)
  ),
  t3 = list(
    data_params = c(1, 3),
    alpha = c(stationary = 0.10, boundary = 0.1508284, explosive = 0.20)
  )
)

core_count <- parallel::detectCores(logical = TRUE)
workers_arg <- NA_integer_
workers_arg <- as.integer(get_arg_value("workers", NA_character_))
workers_env <- suppressWarnings(as.integer(Sys.getenv("MOMENT_POWER_COMPARISON_WORKERS", "")))
workers <- if (!is.na(workers_arg)) {
  workers_arg
} else if (!is.na(workers_env)) {
  workers_env
} else {
  max(1L, floor(core_count * 2 / 3))
}
workers <- max(1L, min(workers, core_count))

message(sprintf(
  "Detected %d logical cores; using %d workers.",
  core_count,
  workers
))

cluster <- parallel::makeCluster(workers)
doParallel::registerDoParallel(cluster)
on.exit({
  try(parallel::stopCluster(cluster), silent = TRUE)
}, add = TRUE)

parallel::clusterExport(
  cluster,
  varlist = c("old_workspace"),
  envir = environment()
)

parallel::clusterEvalQ(cluster, {
  suppressPackageStartupMessages({
    library(gamlss.dist)
    library(stabledist)
  })
  load(old_workspace, envir = .GlobalEnv)
  NULL
})

cells <- do.call(
  rbind,
  lapply(names(innovation_specs), function(innovation) {
    spec <- innovation_specs[[innovation]]
    regimes <- intersect(regime_grid, names(spec$alpha))
    if (length(regimes) == 0L) {
      stop(sprintf("No valid regimes requested for innovation %s.", innovation))
    }
    expand.grid(
      innovation = innovation,
      regime = regimes,
      z0 = z_grid,
      replication = seq_len(replications),
      stringsAsFactors = FALSE
    )
  })
)
cells$cell_id <- seq_len(nrow(cells))

run_replication <- function(task) {
  innovation <- task[["innovation"]]
  regime <- task[["regime"]]
  z0 <- as.numeric(task[["z0"]])
  replication <- as.integer(task[["replication"]])
  cell_id <- as.integer(task[["cell_id"]])

  spec <- innovation_specs[[innovation]]
  alpha0 <- unname(spec$alpha[[regime]])
  data_params <- spec$data_params
  seed <- seed_base + cell_id

  result <- tryCatch({
    set.seed(seed)
    eta <- reta(n_samples + 1000L, params = data_params)

    theta_dgp <- c(0.01, alpha0, 0.9, z0)
    yt <- rGARCH2(n_samples, theta_dgp, eta, m = 1)

    ## Historical TS_NGQMLE_func reads this global optimizer start.
    theta0 <<- c(0.1, 0.1, 0.6, 0.1)

    fit <- TS_NGQMLE_process(
      yt,
      r_power,
      ngqmle_params,
      1,
      method = 1,
      Adaptive = 1
    )

    theta_hat <- fit$TS_NGQMLE1$hat_theta
    lambda_hat <- fit$lambda_hat1
    qmle_params <- fit$qmle_params

    score_test <- Sigma_h(
      m_lag,
      s_power,
      r_power,
      yt,
      theta_hat,
      qmle_params,
      lambda_hat
    )
    li_test <- Sigma_rs(
      m_lag,
      s_power,
      r_power,
      yt,
      theta_hat,
      qmle_params,
      lambda_hat
    )

    score_p <- score_test$pvalue[m_lag]
    li_p <- li_test$pvalue[m_lag]

    data.frame(
      innovation = innovation,
      regime = regime,
      n = n_samples,
      z0 = z0,
      replication = replication,
      alpha0 = alpha0,
      p_score = score_p,
      p_li = li_p,
      reject_score = as.integer(score_p < 0.05),
      reject_li = as.integer(li_p < 0.05),
      status = "ok",
      error = NA_character_,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      innovation = innovation,
      regime = regime,
      n = n_samples,
      z0 = z0,
      replication = replication,
      alpha0 = alpha0,
      p_score = NA_real_,
      p_li = NA_real_,
      reject_score = NA_integer_,
      reject_li = NA_integer_,
      status = "error",
      error = conditionMessage(e),
      stringsAsFactors = FALSE
    )
  })

  result
}

message(sprintf("Running %d Monte Carlo tasks...", nrow(cells)))

raw_results <- foreach(
  i = seq_len(nrow(cells)),
  .combine = rbind,
  .inorder = FALSE,
  .packages = c("gamlss.dist", "stabledist"),
  .export = c(
    "run_replication",
    "innovation_specs",
    "n_samples",
    "m_lag",
    "r_power",
    "s_power",
    "ngqmle_params",
    "seed_base"
  )
) %dopar% {
  run_replication(cells[i, ])
}

write.csv(
  raw_results,
  file.path(out_dir, "raw_results.csv"),
  row.names = FALSE
)

aggregate_cell <- function(df) {
  ok <- df[df$status == "ok", , drop = FALSE]
  data.frame(
    innovation = df$innovation[1],
    regime = df$regime[1],
    n = df$n[1],
    z0 = df$z0[1],
    replications = nrow(df),
    ok = nrow(ok),
    errors = sum(df$status != "ok"),
    rejection_score = mean(ok$reject_score),
    rejection_li = mean(ok$reject_li),
    se_score = sqrt(mean(ok$reject_score) * (1 - mean(ok$reject_score)) / nrow(ok)),
    se_li = sqrt(mean(ok$reject_li) * (1 - mean(ok$reject_li)) / nrow(ok)),
    stringsAsFactors = FALSE
  )
}

split_key <- paste(
  raw_results$innovation,
  raw_results$regime,
  raw_results$n,
  raw_results$z0,
  sep = "|"
)
summary_stats <- do.call(rbind, lapply(split(raw_results, split_key), aggregate_cell))
summary_stats <- summary_stats[order(
  summary_stats$innovation,
  match(summary_stats$regime, c("stationary", "boundary", "explosive")),
  summary_stats$z0
), ]

write.csv(
  summary_stats,
  file.path(out_dir, "summary_stats.csv"),
  row.names = FALSE
)

make_latex_table <- function(summary_stats) {
  z_values <- sort(unique(summary_stats$z0))
  regimes <- intersect(c("stationary", "boundary", "explosive"), unique(summary_stats$regime))
  wide_rows <- list()
  row_idx <- 1L

  for (innovation in unique(summary_stats$innovation)) {
    for (regime in regimes) {
      row <- summary_stats[
        summary_stats$innovation == innovation & summary_stats$regime == regime,
        ,
        drop = FALSE
      ]
      row <- row[order(row$z0), ]
      wide_row <- data.frame(
        Innovation = innovation,
        Baseline = regime,
        stringsAsFactors = FALSE
      )
      for (z_value in z_values) {
        z_label <- gsub("\\.", "p", sprintf("%.1f", z_value))
        wide_row[[paste0("z", z_label, "_Score")]] <-
          row$rejection_score[abs(row$z0 - z_value) < 1e-12]
        wide_row[[paste0("z", z_label, "_Li")]] <-
          row$rejection_li[abs(row$z0 - z_value) < 1e-12]
      }
      wide_rows[[row_idx]] <- wide_row
      row_idx <- row_idx + 1L
    }
  }

  wide <- do.call(rbind, wide_rows)
  fmt <- function(x) sprintf("%.3f", x)
  tabular_spec <- paste0("ll", paste(rep("cc", length(z_values)), collapse = ""))
  header <- paste(
    "Innovation & Baseline",
    paste(sprintf("\\multicolumn{2}{c}{$z_0=%.1f$}", z_values), collapse = " & "),
    sep = " & "
  )
  header <- paste0(gsub("z_0=0.0", "z_0=0", header, fixed = TRUE), " \\\\")
  cmidrules <- paste(
    sprintf("\\cmidrule(lr){%d-%d}", seq(3, by = 2, length.out = length(z_values)), seq(4, by = 2, length.out = length(z_values))),
    collapse = ""
  )
  subheader <- paste0(" & & ", paste(rep("Score & Li", length(z_values)), collapse = " & "), " \\\\")
  lines <- c(
    "\\begin{table}[htbp]",
    "\\centering",
    "\\caption{Comparison with moment-based residual power diagnostics on the same GNGQMLE residuals.}",
    "\\label{tab:moment_power_comparison}",
    sprintf("\\begin{tabular}{%s}", tabular_spec),
    "\\toprule",
    header,
    cmidrules,
    subheader,
    "\\midrule"
  )

  for (i in seq_len(nrow(wide))) {
    values <- unlist(wide[i, setdiff(names(wide), c("Innovation", "Baseline"))], use.names = FALSE)
    lines <- c(
      lines,
      paste0(wide$Innovation[i], " & ", wide$Baseline[i], " & ", paste(fmt(values), collapse = " & "), " \\\\")
    )
  }

  c(
    lines,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  )
}

latex_lines <- make_latex_table(summary_stats)
writeLines(latex_lines, file.path(out_dir, "moment_power_comparison_table.tex"))

message("Finished. Outputs:")
message(file.path(out_dir, "raw_results.csv"))
message(file.path(out_dir, "summary_stats.csv"))
message(file.path(out_dir, "moment_power_comparison_table.tex"))
print(summary_stats)
