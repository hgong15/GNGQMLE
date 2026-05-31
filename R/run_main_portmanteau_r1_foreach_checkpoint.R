## Checkpointed full main-text portmanteau simulation for Q_1(6), model 5.2/5.3.
##
## This version writes each design-cell result immediately to disk.  If a
## long PSOCK run is interrupted, rerun the script and it will skip complete
## cell files, then combine all completed cells at the end.

suppressPackageStartupMessages({
  library(foreach)
  library(doParallel)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg_value <- function(name, default = NA_character_) {
  pattern <- paste0("^--", name, "=")
  matched <- grep(pattern, args, value = TRUE)
  if (length(matched) == 0L) {
    return(default)
  }
  sub(pattern, "", matched[1])
}

old_workspace <- "R/legacy/legacy_workspace.RData"
replications <- as.integer(get_arg_value("replications", "500"))
model_key <- get_arg_value("model", "52")
if (!model_key %in% c("52", "53")) {
  stop("--model must be either 52 or 53")
}
beta0_dgp <- if (model_key == "52") 0.9 else 0
model_name <- if (model_key == "52") "model 5.2" else "model 5.3"
figure_prefix <- if (model_key == "52") "test3_1" else "test3_2"
out_dir <- get_arg_value(
  "outdir",
  sprintf("results/simulation/portmanteau/main_portmanteau_r1_model%s_rep%d", model_key, replications)
)
cell_dir <- file.path(out_dir, "cell_results")
dir.create(cell_dir, recursive = TRUE, showWarnings = FALSE)

n_grid <- c(500L, 1000L)
z_grid <- seq(0, 0.7, by = 0.1)
m_lag <- 6L
r_power <- 1
s_power <- 1
ngqmle_params <- c(1, 7)
seed_base <- 20260419L

innovation_specs <- list(
  norm = list(
    label = "N(0,1)",
    data_params = c(0, 2),
    alpha = c(stationary = 0.06, boundary = 0.1096508, explosive = 0.16)
  ),
  t5 = list(
    label = "t5",
    data_params = c(1, 5),
    alpha = c(stationary = 0.08, boundary = 0.1201453, explosive = 0.18)
  ),
  gg1 = list(
    label = "gg1",
    data_params = c(0, 1),
    alpha = c(stationary = 0.08, boundary = 0.1206941, explosive = 0.18)
  ),
  t3 = list(
    label = "t3",
    data_params = c(1, 3),
    alpha = c(stationary = 0.10, boundary = 0.1508284, explosive = 0.20)
  )
)

format_z <- function(z) {
  gsub("\\.", "p", sprintf("%.1f", z))
}

cell_key <- function(cell) {
  paste(
    cell[["innovation"]],
    cell[["regime"]],
    paste0("n", as.integer(cell[["n"]])),
    paste0("z", format_z(as.numeric(cell[["z0"]]))),
    sep = "_"
  )
}

cell_path <- function(cell) {
  file.path(cell_dir, paste0(cell_key(cell), ".csv"))
}

cell_is_complete <- function(path) {
  if (!file.exists(path)) {
    return(FALSE)
  }
  df <- tryCatch(read.csv(path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(df)) {
    return(FALSE)
  }
  nrow(df) == replications && all(df$status == "ok")
}

cells <- expand.grid(
  innovation = names(innovation_specs),
  regime = c("stationary", "boundary", "explosive"),
  n = n_grid,
  z0 = z_grid,
  stringsAsFactors = FALSE
)
cells$cell_id <- seq_len(nrow(cells))
cells$path <- vapply(seq_len(nrow(cells)), function(i) cell_path(cells[i, ]), character(1))

missing_idx <- which(!vapply(cells$path, cell_is_complete, logical(1)))
missing_cells <- cells[missing_idx, , drop = FALSE]

core_count <- parallel::detectCores(logical = TRUE)
workers_arg <- as.integer(get_arg_value("workers", NA_character_))
workers_env <- suppressWarnings(as.integer(Sys.getenv("MAIN_PORTMANTEAU_WORKERS", "")))
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
message(sprintf(
  "Complete cells: %d / %d. Missing cells: %d.",
  nrow(cells) - nrow(missing_cells),
  nrow(cells),
  nrow(missing_cells)
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

run_one_replication <- function(innovation, regime, n_samples, z0, replication, seed) {
  spec <- innovation_specs[[innovation]]
  alpha0 <- unname(spec$alpha[[regime]])

  for (attempt in seq_len(50L)) {
    result <- tryCatch({
      set.seed(seed + attempt)
      eta <- reta(n_samples + 1000L, params = spec$data_params)
      theta_dgp <- c(0.01, alpha0, beta0_dgp, z0)
      yt <- rGARCH2(n_samples, theta_dgp, eta, m = 1)

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
      NULL
    })

    if (!is.null(result)) {
      return(result)
    }
  }

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
    error = "failed after 50 attempts",
    stringsAsFactors = FALSE
  )
}

run_cell_to_file <- function(cell) {
  output_path <- cell[["path"]]
  if (file.exists(output_path)) {
    existing <- tryCatch(read.csv(output_path, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(existing) && nrow(existing) == replications && all(existing$status == "ok")) {
      return(data.frame(path = output_path, status = "skipped", stringsAsFactors = FALSE))
    }
  }

  innovation <- cell[["innovation"]]
  regime <- cell[["regime"]]
  n_samples <- as.integer(cell[["n"]])
  z0 <- as.numeric(cell[["z0"]])
  cell_id <- as.integer(cell[["cell_id"]])

  rows <- vector("list", replications)
  for (replication in seq_len(replications)) {
    rows[[replication]] <- run_one_replication(
      innovation = innovation,
      regime = regime,
      n_samples = n_samples,
      z0 = z0,
      replication = replication,
      seed = seed_base + cell_id * 100000L + replication * 100L
    )
  }

  result <- do.call(rbind, rows)
  tmp_path <- paste0(output_path, ".tmp_", Sys.getpid())
  write.csv(result, tmp_path, row.names = FALSE)
  if (file.exists(output_path)) {
    unlink(output_path)
  }
  file.rename(tmp_path, output_path)

  data.frame(
    path = output_path,
    status = if (all(result$status == "ok")) "ok" else "has_errors",
    stringsAsFactors = FALSE
  )
}

if (nrow(missing_cells) > 0L) {
  message(sprintf(
    "Running %d missing design cells x %d replications = %d Monte Carlo fits...",
    nrow(missing_cells),
    replications,
    nrow(missing_cells) * replications
  ))

  status <- foreach(
    i = seq_len(nrow(missing_cells)),
    .combine = rbind,
    .inorder = FALSE,
    .packages = c("gamlss.dist", "stabledist"),
    .export = c(
      "run_cell_to_file",
      "run_one_replication",
      "innovation_specs",
      "replications",
      "m_lag",
      "r_power",
      "s_power",
      "ngqmle_params",
      "seed_base",
      "beta0_dgp"
    )
  ) %dopar% {
    run_cell_to_file(missing_cells[i, ])
  }
  write.csv(status, file.path(out_dir, "cell_status_last_run.csv"), row.names = FALSE)
}

complete_flags <- vapply(cells$path, cell_is_complete, logical(1))
if (!all(complete_flags)) {
  incomplete <- cells[!complete_flags, c("innovation", "regime", "n", "z0", "path")]
  write.csv(incomplete, file.path(out_dir, "incomplete_cells.csv"), row.names = FALSE)
  stop(sprintf(
    "%d cells are still incomplete. Rerun this script to continue.",
    nrow(incomplete)
  ))
}

message("All cells complete. Combining results...")

raw_results <- do.call(rbind, lapply(cells$path, read.csv, stringsAsFactors = FALSE))
write.csv(raw_results, file.path(out_dir, "raw_results.csv"), row.names = FALSE)

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
  match(summary_stats$innovation, names(innovation_specs)),
  match(summary_stats$regime, c("stationary", "boundary", "explosive")),
  summary_stats$n,
  summary_stats$z0
), ]
write.csv(summary_stats, file.path(out_dir, "summary_stats.csv"), row.names = FALSE)

plot_panel <- function(panel_regime, file_path) {
  panel <- summary_stats[summary_stats$regime == panel_regime, , drop = FALSE]
  colors <- c(`500` = "blue", `1000` = "red")
  ltys <- c(`500` = 1, `1000` = 2)
  pchs <- c(norm = 0, t5 = 1, gg1 = 2, t3 = 3)

  pdf(file_path, width = 3, height = 3)
  old_par <- par(mfrow = c(1, 1), mar = c(2, 2, 0.5, 0.5))
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)

  plot(
    NA,
    xlim = range(z_grid),
    ylim = c(0, 1),
    xlab = "",
    ylab = "",
    xaxt = "n"
  )
  axis(1, at = z_grid, labels = sprintf("%.1f", z_grid), cex.axis = 0.8)
  axis(2, cex.axis = 0.8)

  for (innovation in names(innovation_specs)) {
    for (n_value in n_grid) {
      curve <- panel[
        panel$innovation == innovation & panel$n == n_value,
        ,
        drop = FALSE
      ]
      curve <- curve[order(curve$z0), ]
      lines(
        curve$z0,
        curve$rejection_score,
        type = "l",
        col = colors[as.character(n_value)],
        lty = ltys[as.character(n_value)]
      )
      lines(
        curve$z0,
        curve$rejection_score,
        type = "p",
        col = colors[as.character(n_value)],
        pch = pchs[innovation]
      )
    }
  }
  abline(h = 0.05)
}

regime_order <- c("stationary", "boundary", "explosive")
for (j in seq_along(regime_order)) {
  plot_panel(
    regime_order[j],
    file.path(out_dir, sprintf("%s_%d.pdf", figure_prefix, j))
  )
}

writeLines(
  c(
    "# Main portmanteau simulation",
    "",
    sprintf("- Date: %s", Sys.time()),
    sprintf("- DGP: %s.", model_name),
    sprintf("- beta0 in the DGP: %.1f.", beta0_dgp),
    sprintf("- Replications per cell: %d", replications),
    "- Estimator/test: GNGQMLE residuals with r=1; score test Q_1(6).",
    "- Innovations: norm, t5, gg1, t3.",
    "- Baselines: stationary, boundary, explosive.",
    "- n: 500, 1000.",
    "- z0: 0, 0.1, ..., 0.7.",
    "",
    "Outputs:",
    "- raw_results.csv",
    "- summary_stats.csv",
    sprintf("- %s_1.pdf, %s_2.pdf, %s_3.pdf", figure_prefix, figure_prefix, figure_prefix)
  ),
  file.path(out_dir, "RESULTS_NOTE.md")
)

message("Finished. Outputs:")
message(file.path(out_dir, "raw_results.csv"))
message(file.path(out_dir, "summary_stats.csv"))
message(file.path(out_dir, sprintf("%s_1.pdf", figure_prefix)))
message(file.path(out_dir, sprintf("%s_2.pdf", figure_prefix)))
message(file.path(out_dir, sprintf("%s_3.pdf", figure_prefix)))
print(summary_stats)
