## Cryptocurrency application diagnostics for the Supplementary Material.
##
## This script computes the same stationarity and score-portmanteau quantities
## as Table 5 in the main text:
##   gamma_hat, Tn, Q_1(6), Q_1(12), and min_{m<=15} p_m.

load_old_application_functions <- function() {
  suppressPackageStartupMessages({
    library(gamlss.dist)
    library(PearsonDS)
  })

  dsgg <<- function(x, beta, nu, lambda, log = FALSE) {
    c0 <- beta / 2^(1 + 1 / beta) / gamma(1 / beta)
    y <- -(abs(x) / lambda / (1 + nu * sign(x)))^beta / 2 -
      log(lambda) + log(c0)
    if (log) y else exp(y)
  }

  source("R/legacy/useful_fun2_subset.R", local = .GlobalEnv)

  source("R/legacy/Eh1h2Er1r2_subset.R", local = .GlobalEnv)
}

load_old_application_functions()

use_reported_crypto_estimates <- TRUE

reported_crypto_fit <- function(series) {
  fits <- list(
    BTC = list(
      theta = c(3.28127161705237e-05, 0.0713967724898285, 0.827487191124844),
      lambda = 1.45475840685351,
      qmle_params = c(0, 0.898077751690437),
      llf = 2835.33850846662,
      aic = -5660.7
    ),
    ETH = list(
      theta = c(1.6e-04, 0.0852, 0.7301),
      lambda = 1.4476,
      qmle_params = c(0, 0.9135),
      llf = 2383.1,
      aic = -4756.2
    ),
    BNB = list(
      theta = c(1.4e-04, 0.0882, 0.7510),
      lambda = 0.9518,
      qmle_params = c(1, 3.3463),
      llf = 2333.9,
      aic = -4657.8
    ),
    TRX = list(
      theta = c(9.1e-05, 0.0960, 0.7845),
      lambda = 0.9286,
      qmle_params = c(1, 3.1120),
      llf = 2230.0,
      aic = -4450.1
    )
  )
  fits[[series]]
}

lambda_estimate_func2_symmetric <- function(error_tilde) {
  L0 <- function(lambda) {
    theta1 <- exp(lambda)
    -mean(dPE(error_tilde / theta1[1], nu = theta1[2], log = TRUE)) +
      log(theta1[1])
  }
  opt0 <- optim(c(0, 0), L0)

  L1 <- function(lambda) {
    theta1 <- exp(lambda)
    -mean(log(dt(error_tilde / theta1[1], df = theta1[2]) / theta1[1]))
  }
  opt1 <- optim(c(0, 1), L1)

  if (opt0$value <= opt1$value) {
    return(c(0, exp(opt0$par)))
  }
  c(1, exp(opt1$par))
}

read_crypto_returns <- function(symbol, start_date, end_date) {
  path <- sprintf("data/raw/%s-usd-max.csv",
                  tolower(symbol))
  if (!file.exists(path)) {
    stop("Cannot find cryptocurrency data file: ", path)
  }
  dat <- read.csv(path, stringsAsFactors = FALSE)
  dat$date <- as.Date(substr(dat$snapped_at, 1, 10))
  dat$price <- as.numeric(dat$price)
  dat <- dat[is.finite(dat$price) & !is.na(dat$date), ]
  dat <- dat[order(dat$date), ]

  y <- diff(log(dat$price))
  date <- dat$date[-1]
  keep <- date >= as.Date(start_date) & date <= as.Date(end_date)
  as.numeric(y[keep])
}

window_stats <- function(y) {
  z <- (y - mean(y)) / sd(y)
  data.frame(
    n = length(y),
    mean = mean(y),
    sd = sd(y),
    skewness = mean(z^3),
    kurtosis = mean(z^4)
  )
}

eta_theta_derivatives <- function(yt, theta) {
  omega <- theta[1]
  alpha <- theta[2]
  beta <- theta[3]
  n <- length(yt)
  ht <- rep(0, n)
  ht[1] <- 0.1
  eta <- rep(0, n)
  psimga_po <- rep(0, n)
  psimga_pa <- rep(0, n)
  psimga_pb <- rep(0, n)

  for (i in 2:n) {
    ht[i] <- omega + alpha * yt[i - 1]^2 + beta * ht[i - 1]
    eta[i] <- yt[i] / sqrt(ht[i])
    psimga_po[i] <- psimga_po[i - 1] + beta^(i - 1)

    betav <- exp(seq(i - 2, 0, by = -1) * log(beta))
    yt2 <- yt[1:(i - 1)]^2
    psimga_pa[i] <- sum(betav * yt2)

    if (i >= 3) {
      betav2 <- exp(seq(i - 3, 0, by = -1) * log(beta))
      psimga_pb[i] <- sum(seq(i - 2, 1, -1) * betav2 *
                            (omega + alpha * yt[1:(i - 2)]^2))
    }
  }

  list(
    hat_eta = eta,
    ht = ht,
    psimga_ptheta = list(
      psimga_po = psimga_po / ht,
      psimga_pa = psimga_pa / ht,
      psimga_pb = psimga_pb / ht
    )
  )
}

kappa_r_local <- function(eta, r) {
  (mean(abs(eta)^(2 * r), na.rm = TRUE) - 1) * 4 / r^2
}

kappa_f_local <- function(eta, lambda, qmle_params) {
  family <- qmle_params[1]
  if (family == 0) {
    return(4 * mean(hgg_1(eta, qmle_params[2], lambda)^2, na.rm = TRUE) /
             lambda^2 /
             mean(hgg_2(eta, qmle_params[2], lambda), na.rm = TRUE)^2)
  }
  if (family == 1) {
    return(4 * mean(ht_1(eta, qmle_params[2], lambda)^2, na.rm = TRUE) /
             lambda^2 /
             mean(ht_2(eta, qmle_params[2], lambda), na.rm = TRUE)^2)
  }
  if (family == 2) {
    return(4 * mean(hPIV_1(eta, qmle_params[2], qmle_params[3], lambda)^2,
                    na.rm = TRUE) / lambda^2 /
             mean(hPIV_2(eta, qmle_params[2], qmle_params[3], lambda),
                  na.rm = TRUE)^2)
  }
  if (family == 3) {
    return(4 * mean(hsgg_1(eta, qmle_params[2], qmle_params[3], lambda)^2,
                    na.rm = TRUE) / lambda^2 /
             mean(hsgg_2(eta, qmle_params[2], qmle_params[3], lambda),
                  na.rm = TRUE)^2)
  }
  stop("Unsupported quasi-likelihood family code: ", family)
}

score_components <- function(eta, qmle_params, lambda) {
  family <- qmle_params[1]
  if (family == 0) {
    return(list(
      h1 = hgg_1(eta, qmle_params[2], lambda),
      h2 = hgg_2(eta, qmle_params[2], lambda)
    ))
  }
  if (family == 1) {
    return(list(
      h1 = ht_1(eta, qmle_params[2], lambda),
      h2 = ht_2(eta, qmle_params[2], lambda)
    ))
  }
  if (family == 2) {
    return(list(
      h1 = hPIV_1(eta, qmle_params[2], qmle_params[3], lambda),
      h2 = hPIV_2(eta, qmle_params[2], qmle_params[3], lambda)
    ))
  }
  if (family == 3) {
    return(list(
      h1 = hsgg_1(eta, qmle_params[2], qmle_params[3], lambda),
      h2 = hsgg_2(eta, qmle_params[2], qmle_params[3], lambda)
    ))
  }
  stop("Unsupported quasi-likelihood family code: ", family)
}

score_portmanteau <- function(m, r, yt, theta, qmle_params, lambda) {
  n <- length(yt)
  result <- eta_theta_derivatives(yt, theta)
  eta <- result$hat_eta
  sc <- score_components(eta, qmle_params, lambda)

  h1_eta <- sc$h1 - mean(sc$h1, na.rm = TRUE)
  s2_h1 <- var(h1_eta, na.rm = TRUE)
  psimga_ptheta <- result$psimga_ptheta

  pmat <- matrix(NA_real_, m, 3)
  rho <- matrix(NA_real_, m, 1)
  for (k in 1:m) {
    pmat[k, ] <- c(
      mean(h1_eta[1:(n - k)] * psimga_ptheta$psimga_po[(k + 1):n],
           na.rm = TRUE),
      mean(h1_eta[1:(n - k)] * psimga_ptheta$psimga_pa[(k + 1):n],
           na.rm = TRUE),
      mean(h1_eta[1:(n - k)] * psimga_ptheta$psimga_pb[(k + 1):n],
           na.rm = TRUE)
    )
    rho[k, 1] <- mean(h1_eta[1:(n - k)] * h1_eta[(k + 1):n],
                      na.rm = TRUE) / sqrt(n - k) / s2_h1
  }

  ka_r <- kappa_r_local(eta, r)
  ka_f <- kappa_f_local(eta, lambda, qmle_params)
  s1 <- mean(psimga_ptheta$psimga_po^2, na.rm = TRUE)
  sigma_23_1 <- matrix(c(
    mean(psimga_ptheta$psimga_po * psimga_ptheta$psimga_pa, na.rm = TRUE),
    mean(psimga_ptheta$psimga_po * psimga_ptheta$psimga_pb, na.rm = TRUE)
  ), nrow = 2)
  sigma_23_23 <- matrix(c(
    mean(psimga_ptheta$psimga_pa^2, na.rm = TRUE),
    mean(psimga_ptheta$psimga_pa * psimga_ptheta$psimga_pb, na.rm = TRUE),
    0,
    mean(psimga_ptheta$psimga_pb^2, na.rm = TRUE)
  ), nrow = 2, byrow = TRUE)
  sigma_23_23[2, 1] <- sigma_23_23[1, 2]

  jmat <- matrix(NA_real_, 3, 3)
  jmat[1, 1] <- s1
  jmat[2:3, 1] <- sigma_23_1
  jmat[1, 2:3] <- t(sigma_23_1)
  jmat[2:3, 2:3] <- sigma_23_23

  b1 <- matrix(c(theta[1], theta[2], 0), nrow = 3)
  jinv <- solve(jmat)
  sigma1 <- ka_f * jinv + (ka_r - ka_f) * (b1 %*% t(b1))
  sigma0 <- diag(rep(1, m)) -
    2 / s2_h1 * pmat %*% jinv %*% t(pmat) +
    1 / ka_f / s2_h1 * pmat %*% sigma1 %*% t(pmat)

  q_stat <- rep(NA_real_, m)
  pvalue <- rep(NA_real_, m)
  for (i in 1:m) {
    q_stat[i] <- n * (n + 2) *
      t(rho[1:i]) %*% solve(sigma0[1:i, 1:i]) %*% rho[1:i]
    pvalue[i] <- 1 - pchisq(q_stat[i], i)
  }
  list(Qm = q_stat, pvalue = pvalue, ka_r = ka_r, ka_f = ka_f)
}

gamma_stat <- function(y, theta) {
  eta <- eta_theta_derivatives(y, theta)$hat_eta
  x1 <- log(theta[2] * eta^2 + theta[3])
  gamma_hat <- mean(x1)
  sigma_u_hat <- mean(x1^2) - gamma_hat^2
  Tn <- sqrt(length(y)) * gamma_hat / sqrt(sigma_u_hat)
  c(gamma_hat = gamma_hat, Tn = Tn)
}

fit_gngqmle_r1 <- function(y) {
  r <- 1
  initial_family <- c(1, 7)
  fit0 <- QMLEstimator_func(y, r = r, method = 1)
  theta0 <<- fit0$QMLEstimator
  residuals0 <- eta_theta(y, theta0)$hat_eta
  lambda0 <- lambda_estimate_func(residuals0, initial_family)
  fit_prelim <- TS_NGQMLE_func(y, lambda0, initial_family, method = 1)
  theta0 <<- fit_prelim$QMLEstimator
  selected <- lambda_estimate_func2_symmetric(residuals0)
  lambda_hat <- selected[2]
  qmle_params <- selected[-2]
  fit1 <- TS_NGQMLE_func(y, lambda_hat, qmle_params, method = 1)

  list(
    theta = fit1$QMLEstimator,
    lambda = lambda_hat,
    qmle_params = qmle_params,
    llf = fit1$log_likelihood,
    aic = fit1$AIC + ifelse(length(qmle_params) == 2, 4, 6)
  )
}

fit_case <- function(series, symbol, start_date, end_date) {
  y <- read_crypto_returns(symbol, start_date, end_date)
  fit <- if (use_reported_crypto_estimates) {
    reported_crypto_fit(series)
  } else {
    fit_gngqmle_r1(y)
  }
  stat <- gamma_stat(y, fit$theta)
  diag <- score_portmanteau(15, 1, y, fit$theta, fit$qmle_params, fit$lambda)

  list(
    summary = cbind(
      series = series,
      start = start_date,
      end = end_date,
      window_stats(y)
    ),
    estimates = data.frame(
      series = series,
      omega = fit$theta[1],
      alpha = fit$theta[2],
      beta = fit$theta[3],
      alpha_plus_beta = fit$theta[2] + fit$theta[3],
      family_code = fit$qmle_params[1],
      lambda = fit$lambda,
      qpar1 = fit$qmle_params[2],
      qpar2 = ifelse(length(fit$qmle_params) >= 3, fit$qmle_params[3], NA),
      llf = fit$llf,
      aic = fit$aic
    ),
    diagnostics_full = data.frame(
      series = series,
      m = 1:15,
      Qm = diag$Qm,
      pvalue = diag$pvalue
    ),
    diagnostics_summary = data.frame(
      series = series,
      gamma_hat = stat[["gamma_hat"]],
      Tn = stat[["Tn"]],
      Q6 = diag$Qm[6],
      p6 = diag$pvalue[6],
      Q12 = diag$Qm[12],
      p12 = diag$pvalue[12],
      min_p_1_15 = min(diag$pvalue),
      below_5pct_1_15 = sum(diag$pvalue < 0.05)
    )
  )
}

cases <- list(
  fit_case("BTC", "btc", "2017-01-01", "2020-12-31"),
  fit_case("ETH", "eth", "2017-01-01", "2020-12-31"),
  fit_case("BNB", "bnb", "2018-01-01", "2021-12-31"),
  fit_case("TRX", "trx", "2018-01-01", "2021-12-31")
)

out_dir <- file.path("results", "application_final")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

summary_stats <- do.call(rbind, lapply(cases, `[[`, "summary"))
estimates <- do.call(rbind, lapply(cases, `[[`, "estimates"))
diagnostics_full <- do.call(rbind, lapply(cases, `[[`, "diagnostics_full"))
diagnostics_summary <- do.call(rbind, lapply(cases, `[[`, "diagnostics_summary"))

write.csv(summary_stats, file.path(out_dir, "crypto_summary_stats.csv"),
          row.names = FALSE)
write.csv(estimates, file.path(out_dir, "crypto_gngqmle_estimates.csv"),
          row.names = FALSE)
write.csv(diagnostics_full, file.path(out_dir, "crypto_diagnostics_full.csv"),
          row.names = FALSE)
write.csv(diagnostics_summary,
          file.path(out_dir, "crypto_diagnostics_summary.csv"),
          row.names = FALSE)

cat("\nSummary statistics:\n")
print(summary_stats, digits = 6)
cat("\nGNGQMLE estimates:\n")
print(estimates, digits = 6)
cat("\nStationarity and diagnostic summary:\n")
print(diagnostics_summary, digits = 6)

cat("\nLaTeX rows for Supplementary Material:\n")
for (i in seq_len(nrow(diagnostics_summary))) {
  row <- diagnostics_summary[i, ]
  cat(sprintf(
    "%s & %.4f & %.2f & %.2f (%.3f) & %.2f (%.3f) & %.3f \\\\\n",
    row$series, row$gamma_hat, row$Tn,
    row$Q6, row$p6, row$Q12, row$p12, row$min_p_1_15
  ))
}
