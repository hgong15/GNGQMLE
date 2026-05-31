## Final exchange-rate application results for the JBES manuscript.

load_old_application_functions <- function() {
  suppressPackageStartupMessages({
    library(gamlss.dist)
    library(numDeriv)
    library(PearsonDS)
  })
  dsgg <<- function(x, beta, nu, lambda, log = FALSE) {
    C <- beta / 2^(1 + 1 / beta) / gamma(1 / beta)
    y <- -(abs(x) / lambda / (1 + nu * sign(x)))^beta / 2 -
      log(lambda) + log(C)
    if (log) y else exp(y)
  }
  source("R/legacy/useful_fun2_subset.R", local = .GlobalEnv)

  source("R/legacy/Eh1h2Er1r2_subset.R", local = .GlobalEnv)
}

load_old_application_functions()

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

kappa_r_local <- function(eta, r) {
  (mean(abs(eta)^(2 * r), na.rm = TRUE) - 1) * 4 / r^2
}

kappa_f_local <- function(eta, lambda, data_params) {
  if (data_params[1] == 0) {
    return(4 * mean(hgg_1(eta, data_params[2], lambda)^2, na.rm = TRUE) /
             lambda^2 / mean(hgg_2(eta, data_params[2], lambda),
                             na.rm = TRUE)^2)
  }
  if (data_params[1] == 1) {
    return(4 * mean(ht_1(eta, data_params[2], lambda)^2, na.rm = TRUE) /
             lambda^2 / mean(ht_2(eta, data_params[2], lambda),
                             na.rm = TRUE)^2)
  }
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
  stop("Unsupported quasi-likelihood family.")
}

Simga1_hat_local <- function(yt, theta, r = 0.5, qmle_params) {
  result <- eta_theta_2(yt, theta)
  eta <- result$hat_eta
  ka_r <- kappa_r_local(eta, r)
  ka_f <- kappa_f_local(eta, qmle_params[1], qmle_params[-1])
  psimga_ptheta <- result$psimga_ptheta
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
  J <- matrix(NA_real_, 3, 3)
  J[1, 1] <- s1
  J[2:3, 1] <- sigma_23_1
  J[1, 2:3] <- t(sigma_23_1)
  J[2:3, 2:3] <- sigma_23_23
  b1 <- matrix(c(theta[1], theta[2], 0), nrow = 3)
  S1 <- ka_f * solve(J) + (ka_r - ka_f) * (b1 %*% t(b1))
  SG <- ka_r * solve(J)
  list(S1 = S1, SG = SG, J = J)
}

log_piv_local <- function(x, phi) {
  lambda <- phi[1]
  m <- phi[2]
  nu <- phi[3]
  if (lambda <= 0 || m <= 0.5) return(rep(-1e100, length(x)))
  out <- PearsonDS::dpearsonIV(
    x, m = m, nu = nu, location = 0, scale = lambda, log = TRUE
  )
  out[!is.finite(out)] <- -1e100
  out
}

log_sgg_local <- function(x, phi) {
  lambda <- phi[1]
  beta <- phi[2]
  nu <- phi[3]
  if (lambda <= 0 || beta <= 0 || abs(nu) >= 1) {
    return(rep(-1e100, length(x)))
  }
  A <- 1 + nu * sign(x)
  if (any(A <= 0)) return(rep(-1e100, length(x)))
  log_c <- log(beta) - (1 + 1 / beta) * log(2) - lgamma(1 / beta)
  out <- -0.5 * (abs(x) / (lambda * A))^beta - log(lambda) + log_c
  out[!is.finite(out)] <- -1e100
  out
}

shape_sigma_local <- function(eta, phi, r, log_density) {
  eta <- eta[is.finite(eta)]
  n <- length(eta)
  score_mat <- numDeriv::jacobian(function(p) log_density(eta, p), phi)
  score_mat <- as.matrix(score_mat)
  n_mat <- crossprod(score_mat) / n
  m_mat <- numDeriv::hessian(function(p) sum(log_density(eta, p)), phi) / n
  b_phi <- (4 / (r * phi[1])) *
    matrix(colMeans((1 - abs(eta)^r) * score_mat), ncol = 1)
  kappa_r <- 4 * (mean(abs(eta)^(2 * r)) - 1) / r^2
  e1 <- matrix(c(1, rep(0, length(phi) - 1)), ncol = 1)
  m_inv <- solve(m_mat)
  sigma_phi <- m_inv %*% n_mat %*% m_inv +
    phi[1]^2 / 4 *
      (kappa_r * e1 %*% t(e1) -
         m_inv %*% b_phi %*% t(e1) -
         e1 %*% t(b_phi) %*% m_inv)
  sigma_phi <- Re(sigma_phi)
  dimnames(sigma_phi) <- list(c("lambda", "nu1", "nu2"),
                              c("lambda", "nu1", "nu2"))
  sigma_phi
}

shape_se_local <- function(yt, theta, r, family, lambda, qpar1, qpar2) {
  eta <- eta_theta_2(yt, theta)$hat_eta
  phi <- c(lambda, qpar1, qpar2)
  log_density <- switch(
    as.character(family),
    "2" = log_piv_local,
    "3" = log_sgg_local,
    stop("Shape standard errors are only implemented for pIV and sgg.")
  )
  sigma_phi <- shape_sigma_local(eta, phi, r, log_density)
  diag_vals <- diag(sigma_phi)
  if (any(diag_vals < 0)) {
    warning("Negative diagonal element in shape covariance matrix.")
  }
  sqrt(ifelse(diag_vals >= 0, diag_vals, NA_real_) / length(eta))
}

TS_NGQMLE_process <- function(yt, r = 2, NGQMLE_params = c(1, 7),
                              Cr = 1, method = 1, method_lambda = 1,
                              Adaptive = 1) {
  result1 <- QMLEstimator_func(yt, r = r, method = method)
  tilde_theta <- result1$QMLEstimator
  theta0 <<- tilde_theta
  result2 <- eta_theta(yt, tilde_theta)
  hat_lambda <- lambda_estimate_func(result2$hat_eta, NGQMLE_params)

  params2 <- NULL
  params3 <- NULL
  if (Adaptive == 1) {
    params2 <- lambda_estimate_func2(result2$hat_eta, method = method_lambda)
    params3 <- params2[-2]
  }

  result3 <- TS_NGQMLE_func(yt, hat_lambda, NGQMLE_params, method = method)
  theta0 <<- result3$QMLEstimator
  result4 <- result3
  if (Adaptive == 1) {
    result4 <- TS_NGQMLE_func(yt, params2[2], params3, method = method)
  }

  list(
    GQMLE = list(
      tilde_theta = c(tilde_theta[1:2] / Cr^2, tilde_theta[3]),
      log_likelihood = result1$log_likelihood,
      AIC = result1$AIC
    ),
    TS_NGQMLE0 = list(
      hat_theta = c(result3$QMLEstimator[1:2] / Cr^2,
                    result3$QMLEstimator[3]),
      log_likelihood = result3$log_likelihood,
      AIC = result3$AIC
    ),
    TS_NGQMLE1 = list(
      hat_theta = c(result4$QMLEstimator[1:2] / Cr^2,
                    result4$QMLEstimator[3]),
      log_likelihood = result4$log_likelihood,
      AIC = result4$AIC
    ),
    lambda_hat0 = hat_lambda,
    lambda_hat1 = params2[2],
    qmle_params = params3
  )
}

read_exchange_returns <- function(symbol, start_date, end_date,
                                  jitter_zero = FALSE) {
  path <- sprintf("data/raw/USD_%s_raw.csv", symbol)[1]
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

window_stats <- function(y) {
  mu <- mean(y)
  sig <- sd(y)
  z <- (y - mu) / sig
  data.frame(
    n = length(y),
    mean = mu,
    sd = sig,
    skewness = mean(z^3),
    kurtosis = mean(z^4)
  )
}

gamma_stat <- function(y, theta) {
  result <- eta_theta_2(y, theta)
  eta <- result$hat_eta
  x1 <- log(theta[2] * eta^2 + theta[3])
  gamma_hat <- mean(x1)
  sigma_u_hat <- mean(x1^2) - gamma_hat^2
  Tn <- sqrt(length(y)) * gamma_hat / sqrt(sigma_u_hat)
  data.frame(
    gamma_hat = gamma_hat,
    Tn = Tn,
    p_reject_nonstationarity = pnorm(Tn),
    p_reject_stationarity = 1 - pnorm(Tn)
  )
}

quantile_loss <- function(y, q, p) {
  mean((p - as.numeric(y < q)) * (y - q), na.rm = TRUE)
}

garch_filter_simple <- function(y, theta) {
  n <- length(y)
  h <- numeric(n)
  eta <- numeric(n)
  h[1] <- mean(head(y^2, min(50, n)), na.rm = TRUE)
  if (!is.finite(h[1]) || h[1] <= 0) h[1] <- var(y, na.rm = TRUE)
  eta[1] <- y[1] / sqrt(h[1])
  for (i in 2:n) {
    h[i] <- theta[1] + theta[2] * y[i - 1]^2 + theta[3] * h[i - 1]
    eta[i] <- y[i] / sqrt(h[i])
  }
  list(h = h, eta = eta)
}

eval_tail <- function(y, theta, split_frac = 0.70) {
  fit <- garch_filter_simple(y, theta)
  n <- length(y)
  split <- floor(split_frac * n)
  train <- 2:split
  test <- (split + 1):n
  eta_train <- fit$eta[train]
  eta_train <- eta_train[is.finite(eta_train)]
  y_test <- y[test]
  h_test <- fit$h[test]
  ok <- is.finite(y_test) & is.finite(h_test) & h_test > 0
  y_test <- y_test[ok]
  h_test <- h_test[ok]

  out <- data.frame(
    qlike = mean(log(h_test) + y_test^2 / h_test, na.rm = TRUE),
    mse = mean((y_test^2 - h_test)^2, na.rm = TRUE),
    mae = mean(abs(y_test^2 - h_test), na.rm = TRUE)
  )
  for (p in c(0.01, 0.05)) {
    q_eta <- as.numeric(quantile(eta_train, p, type = 8, na.rm = TRUE))
    es_eta <- mean(eta_train[eta_train <= q_eta], na.rm = TRUE)
    var <- q_eta * sqrt(h_test)
    es <- es_eta * sqrt(h_test)
    hit <- y_test < var
    fz0 <- -log((p - 1) * es) -
      (y_test - var) * (p - as.numeric(y_test <= var)) / (p * es) +
      var / es
    out[[paste0("VaR", p * 100, "_hit")]] <- mean(hit, na.rm = TRUE)
    out[[paste0("VaR", p * 100, "_qloss")]] <- quantile_loss(y_test, var, p)
    out[[paste0("ES", p * 100, "_fz0_loss")]] <-
      mean(fz0[is.finite(fz0)], na.rm = TRUE)
  }
  out
}

fit_series <- function(label, symbol, start_date, end_date,
                       jitter_zero = FALSE) {
  y <- read_exchange_returns(symbol, start_date, end_date, jitter_zero)
  n <- length(y)
  ngqmle_params <- c(1, 7)

  fit_r1 <- TS_NGQMLE_process(y, r = 1, NGQMLE_params = ngqmle_params,
                              Cr = 1, method = 1, Adaptive = 1)
  fit_r2 <- TS_NGQMLE_process(y, r = 2, NGQMLE_params = ngqmle_params,
                              Cr = 1, method = 1, Adaptive = 1)

  gng_theta <- fit_r1$TS_NGQMLE1$hat_theta
  lq_theta <- fit_r1$GQMLE$tilde_theta
  gq_theta <- fit_r2$GQMLE$tilde_theta

  gng_se <- sqrt(diag(Simga1_hat_local(
    y, gng_theta, 1, c(fit_r1$lambda_hat1, fit_r1$qmle_params)
  )$S1) / n)
  lq_se <- sqrt(diag(Simga1_hat_local(
    y, lq_theta, 1, c(fit_r1$lambda_hat0, ngqmle_params)
  )$SG) / n)
  gq_se <- sqrt(diag(Simga1_hat_local(
    y, gq_theta, 2, c(fit_r2$lambda_hat0, ngqmle_params)
  )$SG) / n)
  shape_se <- shape_se_local(
    y, gng_theta, 1, fit_r1$qmle_params[1], fit_r1$lambda_hat1,
    fit_r1$qmle_params[2], fit_r1$qmle_params[3]
  )

  estimates <- data.frame(
    series = label,
    method = c("GNGQMLE", "GQMLE", "LQMLE"),
    omega = c(gng_theta[1], gq_theta[1], lq_theta[1]),
    omega_se = c(gng_se[1], gq_se[1], lq_se[1]),
    alpha = c(gng_theta[2], gq_theta[2], lq_theta[2]),
    alpha_se = c(gng_se[2], gq_se[2], lq_se[2]),
    beta = c(gng_theta[3], gq_theta[3], lq_theta[3]),
    beta_se = c(gng_se[3], gq_se[3], lq_se[3]),
    alpha_plus_beta = c(gng_theta[2] + gng_theta[3],
                        gq_theta[2] + gq_theta[3],
                        lq_theta[2] + lq_theta[3]),
    lambda = c(fit_r1$lambda_hat1, NA, NA),
    lambda_se = c(shape_se[1], NA, NA),
    qfamily = c(fit_r1$qmle_params[1], NA, NA),
    qpar1 = c(fit_r1$qmle_params[2], NA, NA),
    qpar1_se = c(shape_se[2], NA, NA),
    qpar2 = c(fit_r1$qmle_params[3], NA, NA),
    qpar2_se = c(shape_se[3], NA, NA),
    llf = c(fit_r1$TS_NGQMLE1$log_likelihood,
            fit_r2$GQMLE$log_likelihood,
            fit_r1$GQMLE$log_likelihood),
    aic = c(fit_r1$TS_NGQMLE1$AIC + 6,
            fit_r2$GQMLE$AIC,
            fit_r1$GQMLE$AIC)
  )

  risks <- do.call(rbind, lapply(seq_len(nrow(estimates)), function(i) {
    theta <- as.numeric(estimates[i, c("omega", "alpha", "beta")])
    cbind(series = label, method = estimates$method[i], eval_tail(y, theta))
  }))

  list(
    summary = cbind(series = label, start = start_date, end = end_date,
                    window_stats(y)),
    estimates = estimates,
    stationarity = cbind(series = label, gamma_stat(y, gng_theta)),
    risk = risks,
    y = y,
    gng_theta = gng_theta,
    gng_lambda = fit_r1$lambda_hat1,
    gng_qmle_params = fit_r1$qmle_params
  )
}

cases <- list(
  fit_series("USD/TRY", "TRY", "2021-01-01", "2024-12-31"),
  fit_series("USD/ARS", "ARS", "2016-01-01", "2019-12-31",
             jitter_zero = TRUE)
)

summary <- do.call(rbind, lapply(cases, `[[`, "summary"))
estimates <- do.call(rbind, lapply(cases, `[[`, "estimates"))
stationarity <- do.call(rbind, lapply(cases, `[[`, "stationarity"))
risk <- do.call(rbind, lapply(cases, `[[`, "risk"))

source("R/check_candidate_model_checking.R")
diagnostics <- do.call(rbind, lapply(cases, function(obj) {
  ans <- Sigma_h(15, 1, obj$y, obj$gng_theta, obj$gng_qmle_params,
                 obj$gng_lambda)
  data.frame(
    series = obj$summary$series[1],
    m = 1:15,
    Qm = ans$Qm,
    pvalue = ans$pvalue
  )
}))

diag_summary <- do.call(rbind, lapply(split(diagnostics, diagnostics$series),
                                      function(x) data.frame(
                                        series = x$series[1],
                                        Q6 = x$Qm[x$m == 6],
                                        p6 = x$pvalue[x$m == 6],
                                        Q12 = x$Qm[x$m == 12],
                                        p12 = x$pvalue[x$m == 12],
                                        min_p_1_15 = min(x$pvalue),
                                        below_5pct_1_15 = sum(x$pvalue < 0.05)
                                      )))

out_dir <- file.path("results", "application")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(summary, file.path(out_dir, "exchange_summary_stats.csv"),
          row.names = FALSE)
write.csv(estimates, file.path(out_dir, "exchange_estimates.csv"),
          row.names = FALSE)
write.csv(stationarity, file.path(out_dir, "exchange_stationarity.csv"),
          row.names = FALSE)
write.csv(diagnostics, file.path(out_dir, "exchange_diagnostics_full.csv"),
          row.names = FALSE)
write.csv(diag_summary, file.path(out_dir, "exchange_diagnostics_summary.csv"),
          row.names = FALSE)
write.csv(risk, file.path(out_dir, "exchange_tail_risk.csv"),
          row.names = FALSE)

print(summary)
print(estimates)
print(stationarity)
print(diag_summary)
print(risk)
