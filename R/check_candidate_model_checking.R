## Model-checking diagnostics for the candidate exchange-rate windows.
##
## Computes the old score-based Q_1(m) p-values for:
##   USD/TRY, 2020-01-01--2024-12-31, PIV GNGQMLE fit
##   USD/ARS, 2016-01-01--2019-12-31, SGG GNGQMLE fit

eta_theta_2 <- function(yt, theta) {
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

  for (i in 2:length(yt)) {
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
  list(hat_eta = eta, ht = ht,
       psimga_ptheta = list(psimga_po = psimga_po / ht,
                            psimga_pa = psimga_pa / ht,
                            psimga_pb = psimga_pb / ht))
}

kappa_r <- function(eta, r) {
  (mean(abs(eta)^(2 * r)) - 1) * 4 / r^2
}

hPIV_1 <- function(x, m, nu, lambda) {
  m * (2 * (x / lambda^2 * (x / lambda)) / (1 + (x / lambda)^2)) +
    nu * (x / lambda^2 / (1 + (x / lambda)^2)) - 1 / lambda
}

hPIV_2 <- function(x, m, nu, lambda) {
  -(nu * (x * (2 * lambda) / (lambda^2)^2 / (1 + (x / lambda)^2) -
             x / lambda^2 * (2 * (x / lambda^2 * (x / lambda))) /
             (1 + (x / lambda)^2)^2) +
      m * (2 * (x / lambda^2 * (x / lambda^2) +
                  x * (2 * lambda) / (lambda^2)^2 * (x / lambda)) /
             (1 + (x / lambda)^2) -
             2 * (x / lambda^2 * (x / lambda)) *
             (2 * (x / lambda^2 * (x / lambda))) /
             (1 + (x / lambda)^2)^2) -
      1 / lambda^2)
}

hsgg_1 <- function(x, beta, nu, lambda) {
  (abs(x) / lambda / (1 + nu * sign(x)))^(beta - 1) *
    (beta * (abs(x) / lambda^2 / (1 + nu * sign(x)))) / 2 -
    1 / lambda
}

hsgg_2 <- function(x, beta, nu, lambda) {
  -(((abs(x) / lambda / (1 + nu * sign(x)))^(beta - 1) *
       (beta * (abs(x) * (2 * lambda) / (lambda^2)^2 /
                  (1 + nu * sign(x)))) +
      (abs(x) / lambda / (1 + nu * sign(x)))^((beta - 1) - 1) *
      ((beta - 1) * (abs(x) / lambda^2 / (1 + nu * sign(x)))) *
      (beta * (abs(x) / lambda^2 / (1 + nu * sign(x))))) / 2 -
      1 / lambda^2)
}

kappa_f <- function(eta, lambda, data_params) {
  if (data_params[1] == 2) {
    return(4 * mean(hPIV_1(eta, data_params[2], data_params[3], lambda)^2,
                    na.rm = TRUE) / lambda^2 /
             mean(hPIV_2(eta, data_params[2], data_params[3], lambda),
                  na.rm = TRUE)^2)
  }
  if (data_params[1] == 3) {
    return(4 * mean(hsgg_1(eta, data_params[2], data_params[3], lambda)^2,
                    na.rm = TRUE) / lambda^2 /
             mean(hsgg_2(eta, data_params[2], data_params[3], lambda),
                  na.rm = TRUE)^2)
  }
  stop("Only PIV and SGG are needed here.")
}

Sigma_h <- function(m, r, yt, theta, qmle_params, lambda) {
  n <- length(yt)
  result <- eta_theta_2(yt, theta)
  eta <- result$hat_eta

  if (qmle_params[1] == 2) {
    h1_eta <- hPIV_1(eta, qmle_params[2], qmle_params[3], lambda)
    h2_eta <- hPIV_2(eta, qmle_params[2], qmle_params[3], lambda)
  } else if (qmle_params[1] == 3) {
    h1_eta <- hsgg_1(eta, qmle_params[2], qmle_params[3], lambda)
    h2_eta <- hsgg_2(eta, qmle_params[2], qmle_params[3], lambda)
  } else {
    stop("Only PIV and SGG are needed here.")
  }

  h1_eta <- h1_eta - mean(h1_eta, na.rm = TRUE)
  s2_h1 <- var(h1_eta, na.rm = TRUE)
  psimga_ptheta <- result$psimga_ptheta

  P <- matrix(NA_real_, m, 3)
  rho <- matrix(NA_real_, m, 1)
  for (k in 1:m) {
    P[k, ] <- c(
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

  ka_r <- kappa_r(eta, r)
  ka_f <- kappa_f(eta, lambda, qmle_params)
  s1 <- mean(psimga_ptheta$psimga_po^2)
  sigma_11 <- matrix(s1, 1)
  sigma_23_1 <- matrix(c(mean(psimga_ptheta$psimga_po *
                                psimga_ptheta$psimga_pa, na.rm = TRUE),
                         mean(psimga_ptheta$psimga_po *
                                psimga_ptheta$psimga_pb, na.rm = TRUE)),
                       nrow = 2)
  sigma_23_23 <- matrix(c(
    mean(psimga_ptheta$psimga_pa^2, na.rm = TRUE),
    mean(psimga_ptheta$psimga_pa * psimga_ptheta$psimga_pb, na.rm = TRUE),
    0,
    mean(psimga_ptheta$psimga_pb^2, na.rm = TRUE)
  ), nrow = 2, byrow = TRUE)
  sigma_23_23[2, 1] <- sigma_23_23[1, 2]
  J <- matrix(NA_real_, 3, 3)
  J[1, 1] <- s1
  J[2:3, 1] <- sigma_23_1
  J[1, 2:3] <- t(sigma_23_1)
  J[2:3, 2:3] <- sigma_23_23

  b1 <- matrix(c(theta[1], theta[2], 0), nrow = 3)
  J_inv <- solve(J)
  Sigma1 <- ka_f * J_inv + (ka_r - ka_f) * (b1 %*% t(b1))
  Sigma0 <- diag(rep(1, m)) - 2 / s2_h1 * P %*% J_inv %*% t(P) +
    1 / ka_f / s2_h1 * P %*% Sigma1 %*% t(P)

  Qm <- rep(NA_real_, m)
  pvalue <- rep(NA_real_, m)
  for (i in 1:m) {
    Sigma0sub <- Sigma0[1:i, 1:i]
    rhosub <- rho[1:i]
    QMsub <- n * (n + 2) * t(rhosub) %*% solve(Sigma0sub) %*% rhosub
    Qm[i] <- QMsub[1, 1]
    pvalue[i] <- 1 - pchisq(Qm[i], i)
  }
  list(Qm = Qm, pvalue = pvalue, ka_r = ka_r, ka_f = ka_f)
}

read_loss_returns <- function(symbol, start_date, end_date, jitter_zero = FALSE) {
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

if (sys.nframe() == 0) {
  cases <- list(
    TRY_2020_2024 = list(
      y = read_loss_returns("TRY", "2020-01-01", "2024-12-31"),
      theta = c(7.38132e-07, 0.395592, 0.486740),
      lambda = 1.2647537,
      qmle_params = c(2, 1.7628206, 0.4293863)
    ),
    ARS_2016_2019 = list(
      y = read_loss_returns("ARS", "2016-01-01", "2019-12-31",
                            jitter_zero = TRUE),
      theta = c(3.61259685918918e-06, 0.230756418093114, 0.549960691362869),
      lambda = 0.221794885592111,
      qmle_params = c(3, 0.741457587487163, -0.107279129829082)
    )
  )

  rows <- list()
  for (nm in names(cases)) {
    obj <- cases[[nm]]
    ans <- Sigma_h(15, 1, obj$y, obj$theta, obj$qmle_params, obj$lambda)
    rows[[length(rows) + 1]] <- data.frame(
      series_window = nm,
      m = 1:15,
      Qm = ans$Qm,
      pvalue = ans$pvalue,
      ka_r = ans$ka_r,
      ka_f = ans$ka_f
    )
  }
  res <- do.call(rbind, rows)

  out_dir <- file.path("results", "application")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(res, file.path(out_dir, "candidate_model_checking_pvalues.csv"),
            row.names = FALSE)

  summary <- aggregate(pvalue ~ series_window, res,
                       function(x) c(min = min(x), below_5pct = sum(x < 0.05)))
  print(res)
  print(summary)
}
