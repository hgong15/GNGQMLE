## Combine n = 500 and n = 1000 moment-based residual power diagnostic comparison summaries into one compact table.

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

input_files <- c(
  sprintf("results/simulation/moment_power_comparison/moment_power_comparison_foreach_r1_n500_rep%d/summary_stats.csv", replications),
  sprintf("results/simulation/moment_power_comparison/moment_power_comparison_foreach_r1_n1000_rep%d/summary_stats.csv", replications)
)
out_dir <- sprintf("results/simulation/moment_power_comparison/moment_power_comparison_foreach_r1_combined_rep%d", replications)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

summary_stats <- do.call(rbind, lapply(input_files, read.csv, stringsAsFactors = FALSE))
summary_stats <- summary_stats[order(
  summary_stats$innovation,
  match(summary_stats$regime, c("stationary", "boundary", "explosive")),
  summary_stats$n,
  summary_stats$z0
), ]

write.csv(
  summary_stats,
  file.path(out_dir, "summary_stats_combined.csv"),
  row.names = FALSE
)

wide_rows <- list()
row_idx <- 1L

for (innovation in unique(summary_stats$innovation)) {
  for (regime in c("stationary", "boundary", "explosive")) {
    for (n_value in sort(unique(summary_stats$n))) {
      row <- summary_stats[
        summary_stats$innovation == innovation &
          summary_stats$regime == regime &
          summary_stats$n == n_value,
        ,
        drop = FALSE
      ]
      row <- row[order(row$z0), ]
      wide_rows[[row_idx]] <- data.frame(
        Innovation = innovation,
        Baseline = regime,
        n = n_value,
        Size_Score = row$rejection_score[row$z0 == 0],
        Size_Li = row$rejection_li[row$z0 == 0],
        Power03_Score = row$rejection_score[row$z0 == 0.3],
        Power03_Li = row$rejection_li[row$z0 == 0.3],
        Power06_Score = row$rejection_score[row$z0 == 0.6],
        Power06_Li = row$rejection_li[row$z0 == 0.6]
      )
      row_idx <- row_idx + 1L
    }
  }
}

wide <- do.call(rbind, wide_rows)
write.csv(
  wide,
  file.path(out_dir, "moment_power_comparison_wide.csv"),
  row.names = FALSE
)

fmt <- function(x) sprintf("%.3f", x)
latex_lines <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Comparison with moment-based residual power diagnostics on the same GNGQMLE residuals.}",
  "\\label{tab:moment_power_comparison}",
  "\\begin{tabular}{lllcccccc}",
  "\\toprule",
  "Innovation & Baseline & $n$ & \\multicolumn{2}{c}{$z_0=0$} & \\multicolumn{2}{c}{$z_0=0.3$} & \\multicolumn{2}{c}{$z_0=0.6$} \\\\",
  "\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}\\cmidrule(lr){8-9}",
  " & & & Score & Li & Score & Li & Score & Li \\\\",
  "\\midrule"
)

for (i in seq_len(nrow(wide))) {
  latex_lines <- c(
    latex_lines,
    sprintf(
      "%s & %s & %d & %s & %s & %s & %s & %s & %s \\\\",
      wide$Innovation[i],
      wide$Baseline[i],
      wide$n[i],
      fmt(wide$Size_Score[i]),
      fmt(wide$Size_Li[i]),
      fmt(wide$Power03_Score[i]),
      fmt(wide$Power03_Li[i]),
      fmt(wide$Power06_Score[i]),
      fmt(wide$Power06_Li[i])
    )
  )
}

latex_lines <- c(
  latex_lines,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)

writeLines(latex_lines, file.path(out_dir, "moment_power_comparison_table_combined.tex"))
print(wide)
