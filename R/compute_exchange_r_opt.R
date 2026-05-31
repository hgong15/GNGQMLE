read_exchange_returns <- function(symbol, start_date, end_date,
                                  jitter_zero = FALSE) {
  path <- sprintf("data/raw/USD_%s_raw.csv", symbol)
  dat <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  names(dat)[1:2] <- c("date", "close")
  dat$date <- as.Date(dat$date, tryFormats = c("%Y-%m-%d", "%Y/%m/%d"))
  dat$close <- as.numeric(dat$close)
  dat <- dat[is.finite(dat$close) & !is.na(dat$date), ]
  dat <- dat[order(dat$date), ]
  y <- -diff(log(dat$close))
  date <- dat$date[-1]
  keep <- date >= as.Date(start_date) & date <= as.Date(end_date)
  y <- as.numeric(y[keep])
  if (jitter_zero && any(y == 0)) {
    set.seed(1)
    y[y == 0] <- rnorm(sum(y == 0), 0, 0.001 * sd(y))
  }
  y
}

eta_filter <- function(y, theta) {
  n <- length(y)
  h <- numeric(n)
  eta <- numeric(n)
  h[1] <- 0.1
  eta[1] <- y[1] / sqrt(h[1])
  for (i in 2:n) {
    h[i] <- theta[1] + theta[2] * y[i - 1]^2 + theta[3] * h[i - 1]
    eta[i] <- y[i] / sqrt(h[i])
  }
  eta[is.finite(eta)]
}

kappa_hat <- function(eta, r) {
  4 / r^2 * (mean(abs(eta)^(2 * r)) / mean(abs(eta)^r)^2 - 1)
}

r_opt_hat <- function(eta, lower = 0.3, upper = 2) {
  opt <- optimize(function(r) kappa_hat(eta, r),
                  interval = c(lower, upper), tol = 1e-8)
  data.frame(
    r_opt = opt$minimum,
    kappa_opt = opt$objective,
    kappa_0_5 = kappa_hat(eta, 0.5),
    kappa_1 = kappa_hat(eta, 1),
    kappa_2 = kappa_hat(eta, 2)
  )
}

est <- read.csv("results/application/exchange_estimates.csv")
configs <- data.frame(
  series = c("USD/TRY", "USD/ARS"),
  symbol = c("TRY", "ARS"),
  start = c("2021-01-01", "2016-01-01"),
  end = c("2024-12-31", "2019-12-31"),
  jitter = c(FALSE, TRUE)
)

results <- do.call(rbind, lapply(seq_len(nrow(configs)), function(i) {
  cfg <- configs[i, ]
  y <- read_exchange_returns(cfg$symbol, cfg$start, cfg$end, cfg$jitter)
  row <- est[est$series == cfg$series & est$method == "LQMLE", ]
  theta <- as.numeric(row[c("omega", "alpha", "beta")])
  eta <- eta_filter(y, theta)
  cbind(series = cfg$series, n = length(y), r_opt_hat(eta))
}))

out_dir <- file.path("results", "application")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(results, file.path(out_dir, "exchange_r_opt.csv"),
          row.names = FALSE)
print(results, digits = 6)
