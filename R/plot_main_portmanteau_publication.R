## Publication-style plots for the main portmanteau simulation.
##
## Keeps the paper layout unchanged: three separate PDFs for the three
## baseline calibration groups.

args <- commandArgs(trailingOnly = TRUE)
get_arg_value <- function(name, default = NA_character_) {
  pattern <- paste0("^--", name, "=")
  matched <- grep(pattern, args, value = TRUE)
  if (length(matched) == 0L) {
    return(default)
  }
  sub(pattern, "", matched[1])
}

result_dir <- get_arg_value("result_dir", "results/simulation/portmanteau/main_portmanteau_r1_model52_rep500")
figure_dir <- get_arg_value("figure_dir", "figures")
figure_prefix <- get_arg_value("figure_prefix", "test3_1")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

summary_stats <- read.csv(file.path(result_dir, "summary_stats.csv"), stringsAsFactors = FALSE)

innovation_order <- c("norm", "t5", "gg1", "t3")
innovation_labels <- c(
  norm = "N(0,1)",
  t5 = expression(t[5]),
  gg1 = expression(gg[1]),
  t3 = expression(t[3])
)

regime_order <- c("stationary", "boundary", "explosive")
regime_titles <- list(
  stationary = "Group A",
  boundary = "Group B",
  explosive = "Group C"
)

colors <- c(
  norm = "#1F5A93",
  t5 = "#C44E52",
  gg1 = "#2E8B57",
  t3 = "#7A5195"
)
pchs <- c(norm = 15, t5 = 16, gg1 = 17, t3 = 3)
ltys <- c(`500` = 1, `1000` = 2)
line_widths <- c(`500` = 1.9, `1000` = 1.9)
n_grid <- c(500, 1000)
z_grid <- sort(unique(summary_stats$z0))

draw_panel <- function(regime, file_path) {
  panel <- summary_stats[summary_stats$regime == regime, , drop = FALSE]

  pdf(file_path, width = 2.95, height = 3.15, useDingbats = FALSE)
  old_par <- par(
    mar = c(2.72, 2.38, 2.95, 0.42),
    mgp = c(1.62, 0.48, 0),
    tcl = -0.24,
    cex.axis = 1.10,
    cex.lab = 1.18,
    family = "serif"
  )
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)

  plot(
    NA,
    xlim = range(z_grid),
    ylim = c(0, 1),
    xlab = expression(z[0]),
    ylab = "",
    xaxt = "n",
    yaxt = "n",
    bty = "n"
  )

  axis(1, at = z_grid, labels = sprintf("%.1f", z_grid))
  axis(2, at = seq(0, 1, by = 0.2), labels = sprintf("%.1f", seq(0, 1, by = 0.2)), las = 1)

  abline(h = seq(0, 1, by = 0.2), col = "#E2E2E2", lwd = 0.65)
  abline(v = z_grid, col = "#EEEEEE", lwd = 0.5)
  abline(h = 0.05, col = "#333333", lwd = 1.1, lty = 3)
  box(col = "#555555", lwd = 0.9)

  for (innovation in innovation_order) {
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
        col = colors[innovation],
        lty = ltys[as.character(n_value)],
        lwd = line_widths[as.character(n_value)]
      )
      points(
        curve$z0,
        curve$rejection_score,
        col = colors[innovation],
        pch = pchs[innovation],
        cex = if (innovation == "t3") 1.08 else 0.92,
        lwd = 1.45
      )
    }
  }

  title(main = regime_titles[[regime]], line = 1.55, cex.main = 1.70, font.main = 1)
}

for (j in seq_along(regime_order)) {
  draw_panel(
    regime_order[j],
    file.path(figure_dir, sprintf("%s_%d.pdf", figure_prefix, j))
  )
}

## Also save a standalone legend PDF in case it is useful for checking or
## future layout changes.  The current paper caption still explains mapping.
pdf(file.path(result_dir, sprintf("%s_legend.pdf", figure_prefix)), width = 4.5, height = 1.2, useDingbats = FALSE)
old_par <- par(mar = c(0, 0, 0, 0), family = "serif")
plot.new()
legend(
  "center",
  legend = c("N(0,1)", "t5", "gg1", "t3", "n=500", "n=1000"),
  col = c(colors, "#333333", "#333333"),
  pch = c(pchs, NA, NA),
  lty = c(rep(NA, 4), 1, 2),
  lwd = c(rep(NA, 4), 1.25, 1.25),
  horiz = TRUE,
  bty = "n",
  cex = 0.85,
  x.intersp = 0.75,
  y.intersp = 0.8
)
par(old_par)
dev.off()

message("Updated publication-style PDFs:")
message(file.path(figure_dir, sprintf("%s_1.pdf", figure_prefix)))
message(file.path(figure_dir, sprintf("%s_2.pdf", figure_prefix)))
message(file.path(figure_dir, sprintf("%s_3.pdf", figure_prefix)))
