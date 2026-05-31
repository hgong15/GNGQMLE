## Create compact wide summaries for the main portmanteau simulation.

args <- commandArgs(trailingOnly = TRUE)
get_arg_value <- function(name, default = NA_character_) {
  pattern <- paste0("^--", name, "=")
  matched <- grep(pattern, args, value = TRUE)
  if (length(matched) == 0L) {
    return(default)
  }
  sub(pattern, "", matched[1])
}

out_dir <- get_arg_value("outdir", "results/simulation/portmanteau/main_portmanteau_r1_model52_rep500")
summary_path <- file.path(out_dir, "summary_stats.csv")
summary_stats <- read.csv(summary_path, stringsAsFactors = FALSE)

summary_stats <- summary_stats[order(
  match(summary_stats$innovation, c("norm", "t5", "gg1", "t3")),
  match(summary_stats$regime, c("stationary", "boundary", "explosive")),
  summary_stats$n,
  summary_stats$z0
), ]

make_wide <- function(value_col) {
  rows <- list()
  idx <- 1L
  for (innovation in c("norm", "t5", "gg1", "t3")) {
    for (regime in c("stationary", "boundary", "explosive")) {
      for (n_value in c(500, 1000)) {
        row <- summary_stats[
          summary_stats$innovation == innovation &
            summary_stats$regime == regime &
            summary_stats$n == n_value,
          ,
          drop = FALSE
        ]
        row <- row[order(row$z0), ]
        values <- setNames(
          as.list(row[[value_col]]),
          paste0("z", gsub("\\.", "p", sprintf("%.1f", row$z0)))
        )
        rows[[idx]] <- data.frame(
          innovation = innovation,
          regime = regime,
          n = n_value,
          values,
          check.names = FALSE
        )
        idx <- idx + 1L
      }
    }
  }
  do.call(rbind, rows)
}

score_wide <- make_wide("rejection_score")
moment_power_wide <- make_wide("rejection_li")

write.csv(score_wide, file.path(out_dir, "summary_score_wide.csv"), row.names = FALSE)
write.csv(moment_power_wide, file.path(out_dir, "summary_moment_power_wide.csv"), row.names = FALSE)

key_points <- summary_stats[
  summary_stats$z0 %in% c(0, 0.3, 0.6, 0.7),
  c("innovation", "regime", "n", "z0", "rejection_score", "rejection_li", "ok", "errors")
]
write.csv(key_points, file.path(out_dir, "summary_key_points.csv"), row.names = FALSE)

cat("Score rejection frequencies, wide format:\n")
print(score_wide)
cat("\nKey points:\n")
print(key_points)
