## Publication-style application diagnostic plots.
##
## The p-values are recovered from the vector coordinates in the previous
## application PDFs and saved in results/application_model_checking.

data_path <- "results/application/estimated_pvalues_from_pdf.csv"
figure_dir <- "figures"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

pvalues <- read.csv(data_path, stringsAsFactors = FALSE)

draw_panel <- function(series_name, file_path) {
  panel <- pvalues[pvalues$series == series_name, , drop = FALSE]
  panel <- panel[order(panel$m), ]

  pdf(file_path, width = 2.95, height = 2.55, useDingbats = FALSE)
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
    xlim = c(1, 15),
    ylim = c(0, 1),
    xlab = expression(italic(m)),
    ylab = "",
    xaxt = "n",
    yaxt = "n",
    bty = "n"
  )

  x_ticks <- seq(2, 14, by = 2)
  y_ticks <- seq(0, 1, by = 0.2)
  axis(1, at = x_ticks, labels = x_ticks)
  axis(2, at = y_ticks, labels = sprintf("%.1f", y_ticks), las = 1)

  abline(h = y_ticks, col = "#E2E2E2", lwd = 0.65)
  abline(v = x_ticks, col = "#EEEEEE", lwd = 0.5)
  abline(h = 0.05, col = "#333333", lwd = 1.1, lty = 3)
  box(col = "#555555", lwd = 0.9)

  lines(
    panel$m,
    panel$pvalue,
    col = "#1F5A93",
    lty = 1,
    lwd = 1.9
  )
  points(
    panel$m,
    panel$pvalue,
    col = "#1F5A93",
    pch = 15,
    cex = 0.92,
    lwd = 1.45
  )

  title(main = series_name, line = 1.15, cex.main = 1.22, font.main = 1)
}

draw_panel("USD/TRY", file.path(figure_dir, "application_TRY.pdf"))
draw_panel("USD/ARS", file.path(figure_dir, "application_ARS.pdf"))

message("Updated application diagnostic PDFs:")
message(file.path(figure_dir, "application_TRY.pdf"))
message(file.path(figure_dir, "application_ARS.pdf"))
