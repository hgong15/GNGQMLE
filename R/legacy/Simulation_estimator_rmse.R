# =========================================================================
# GNGQMLE Simulation Study (Double-Adaptive & Parallelized)
# =========================================================================
library(doParallel)
library(foreach)
library(PearsonDS) # 濡傛灉浣犱娇鐢ㄤ簡澶栭儴鍖呰纭繚鍔犺浇

# =========================================================================
# 0. 鍏ㄥ眬鎺у埗寮€鍏?(TEST MODE)
# =========================================================================
TEST_MODE <- FALSE
#TEST_MODE <- TRUE
if (TEST_MODE) {
  cat("\n======================================================\n")
  cat("!!! 璀﹀憡: 褰撳墠澶勪簬 TEST_MODE (娴嬭瘯妯″紡) !!!\n")
  cat("姝ｅ湪鍗曠嫭娴嬭瘯鏂板鐨?4 涓噸灏惧绉板垎甯? t7, t9, gg0.8, gg0.4\n")
  cat("======================================================\n\n")
  replications <- 200
  n_values <- c(500, 1000,2000)
  regimes <- c("stationary", "boundary", "explosive")
  dgps <- c( "N(0,1)") 
} else {
  replications <- 1000
  n_values <- c(500, 1000)
  regimes <- c("stationary", "boundary", "explosive")
  # 缁堟瀬澶ф弧璐樀瀹?(鍚墍鏈夊熀鍑嗐€佹柊澧炲帤灏句笌闈炲绉版潃鎵?
 # dgps <- c("N(0,1)", "t9", "t7", "t5", "t3", 
  #          "gg1", "gg0.8", "gg0.4", 
   #         "MixN", "MixT", "pIV1", "sgg1") 
  
  dgps <- c("N(0,1)", "t7", "t5","t3", "gg1", "MixN", "MixT", "pIV1", "sgg1")
}

# =========================================================================
# 1. 杞藉叆搴曞眰鍒嗗竷鍑芥暟 (Distributions)
# =========================================================================
rsgg <- function(n, beta, nu, lambda = 1) {
  sign_z <- sample(c(-1, 1), n, replace = TRUE)
  U <- rgamma(n, shape = 1/beta, scale = 2)
  Z <- sign_z * U^(1/beta)
  X <- ifelse(Z >= 0, (1 - nu) * lambda * Z, (1 + nu) * lambda * Z)
  return(X)
}
dsgg <- function(x, beta, nu, lambda = 1, log = FALSE) {
  C = beta / 2^(1 + 1/beta) / gamma(1/beta)
  y = -(abs(x) / lambda / (1 + nu * sign(x)))^beta / 2 - log(lambda) + log(C)
  if(log) return(y) else return(exp(y))
}

dgg <- function(x, nu, sigma = 1) {
  (nu / (2 * sigma * gamma(1/nu))) * exp(-abs(x/sigma)^nu)
}

# --- 娣峰悎鍒嗗竷 ---
rmixN <- function(n) {
  component <- rbinom(n, 1, 0.1) 
  rnorm(n, mean = 0, sd = ifelse(component == 1, 5, 1))
}
dmixN <- function(x) {
  0.9 * dnorm(x, mean = 0, sd = 1) + 0.1 * dnorm(x, mean = 0, sd = 5)
}
rmixT <- function(n) {
  component <- rbinom(n, 1, 0.1) 
  z <- rt(n, df = 5)
  z * ifelse(component == 1, 5, 1)
}
dmixT <- function(x) {
  0.9 * dt(x, df = 5) + 0.1 * dt(x / 5, df = 5) / 5
}

if(!exists("rpearsonIV")) {
  rpearsonIV <- function(n, ...) rnorm(n)
  dpearsonIV <- function(x, ...) dnorm(x)
}
if(!exists("rPE")) {
  rPE <- function(n, mu, sigma, nu) rnorm(n) 
}

# =========================================================================
# 2. 涓ユ牸鐨勭悊璁烘暟鍊肩Н鍒嗚绠楀櫒 (鍙栦唬钂欑壒鍗℃礇)
# =========================================================================
get_rescale_factors_exact <- function(dgp) {
  d_raw <- function(x) {
    switch(dgp,
           "N(0,1)" = dnorm(x),
           "t9"     = dt(x, 9),
           "t7"     = dt(x, 7),
           "t5"     = dt(x, 5),
           "t3"     = dt(x, 3),
           "gg1"    = dgg(x, nu=1),
           "gg0.8"  = dgg(x, nu=0.8),
           "gg0.4"  = dgg(x, nu=0.4),
           "MixN"   = dmixN(x),
           "MixT"   = dmixT(x),
           "pIV1"   = dpearsonIV(x, m=4, nu=1.2, location=0, scale=1),
           "pIV2"   = dpearsonIV(x, m=3, nu=-0.9, location=0, scale=1),
           "sgg1"   = dsgg(x, beta=0.5, nu=0.5, lambda=1),
           "sgg2"   = dsgg(x, beta=0.4, nu=-0.6, lambda=1))
  }
  d_raw_v <- Vectorize(d_raw)
  
  var_func <- function(x) x^2 * d_raw_v(x)
  raw_m2 <- integrate(var_func, -Inf, Inf, rel.tol=1e-8, stop.on.error=FALSE)$value
  rms <- sqrt(raw_m2)
  
  d_scaled <- function(x) { rms * d_raw_v(x * rms) }
  
  c1_func <- function(x) abs(x) * d_scaled(x)
  C1 <- (integrate(c1_func, -Inf, Inf, rel.tol=1e-8)$value)^2
  
  c2_func <- function(x) x^2 * d_scaled(x)
  C2 <- integrate(c2_func, -Inf, Inf, rel.tol=1e-8)$value
  
  prob_inside <- function(c) {
    integrate(d_scaled, -sqrt(c), sqrt(c), rel.tol=1e-8)$value - 0.5
  }
  C3 <- uniroot(prob_inside, lower=0.0001, upper=10, tol=1e-8)$root
  
  return(list(C1 = C1, C2 = C2, C3 = C3, RMS = rms))
}

# =========================================================================
# 3. 鏁版嵁鐢熸垚鍣?
# =========================================================================
generate_eta <- function(dgp, n, rms) {
  eta_raw <- switch(dgp,
                    "N(0,1)" = rnorm(n),
                    "t9"     = rt(n, df=9),
                    "t7"     = rt(n, df=7),
                    "t5"     = rt(n, df=5),
                    "t3"     = rt(n, df=3),
                    "gg1"    = rPE(n, mu=0, sigma=1, nu=1),
                    "gg0.8"  = rPE(n, mu=0, sigma=1, nu=0.8),
                    "gg0.4"  = rPE(n, mu=0, sigma=1, nu=0.4),
                    "MixN"   = rmixN(n),
                    "MixT"   = rmixT(n),
                    "pIV1"   = rpearsonIV(n, m=4, nu=1.2, location=0, scale=1),
                    "pIV2"   = rpearsonIV(n, m=3, nu=-0.9, location=0, scale=1),
                    "sgg1"   = rsgg(n, beta=0.5, nu=0.5, lambda=1),
                    "sgg2"   = rsgg(n, beta=0.4, nu=-0.6, lambda=1)
  )
  return(eta_raw / rms)
}

# =========================================================================
# 4. 鏍稿績浼拌鍑芥暟 
# =========================================================================
garch_filter_log <- function(w, a, b, y) {
  n <- length(y)
  log_h <- numeric(n)
  
  log_w <- log(w)
  log_a <- log(a)
  log_b <- log(b)
  log_y2 <- 2 * log(abs(y) + 1e-150) 
  
  log_h[1] <- log(max(1e-150, w + a * y[1]^2))
  
  for(t in 2:n) {
    l_term2 <- log_a + log_y2[t-1]
    l_term3 <- log_b + log_h[t-1]
    m <- max(log_w, l_term2, l_term3)
    log_h[t] <- m + log(exp(log_w - m) + exp(l_term2 - m) + exp(l_term3 - m))
  }
  return(log_h) 
}

soft_exp <- function(x, limit = 150) {
  ifelse(x < limit, exp(x), exp(limit) + exp(limit) * (x - limit))
}

est_LQMLE <- function(y) {
  log_y2 <- 2 * log(abs(y) + 1e-150)
  obj <- function(p) {
    w <- p[1]; a <- p[2]; b <- p[3]
    if (w <= 0 || a < 0 || b < 0 || b >= 0.999) return(1e9)
    log_h <- garch_filter_log(w, a, b, y)
    u <- log_y2 - log_h
    u <- pmin(u, 100) 
    mean(0.5 * log_h + exp(0.5 * u))
  }
  # 宸蹭慨鏀? alpha 涓婄晫浠?5 鏀逛负 1
  opt <- try(nlminb(c(0.1, 0.1, 0.7), obj, lower=c(1e-6, 1e-6, 1e-6), upper=c(Inf, 1, 0.999)), silent=TRUE)
  if(inherits(opt, "try-error")) return(c(NA, NA, NA))
  return(opt$par)
}

est_GQMLE <- function(y) {
  log_y2 <- 2 * log(abs(y) + 1e-150)
  obj <- function(p) {
    w <- p[1]; a <- p[2]; b <- p[3]
    if (w <= 0 || a < 0 || b < 0 || b >= 0.999) return(1e9)
    log_h <- garch_filter_log(w, a, b, y)
    u <- log_y2 - log_h
    u <- pmin(u, 100)
    mean(0.5 * log_h + 0.5 * exp(u))
  }
  # 宸蹭慨鏀? alpha 涓婄晫浠?5 鏀逛负 1
  opt <- try(nlminb(c(0.1, 0.1, 0.7), obj, lower=c(1e-6, 1e-6, 1e-6), upper=c(Inf, 1, 0.999)), silent=TRUE)
  if(inherits(opt, "try-error")) return(c(NA, NA, NA))
  return(opt$par)
}

est_LADE <- function(y) {
  log_y2 <- 2 * log(abs(y) + 1e-150)
  obj <- function(p) {
    w <- p[1]; a <- p[2]; b <- p[3]
    if (w <= 0 || a < 0 || b < 0 || b >= 0.999) return(1e9)
    log_h <- garch_filter_log(w, a, b, y)
    u <- log_y2 - log_h
    mean(abs(u))
  }
  # 宸蹭慨鏀? alpha 涓婄晫浠?5 鏀逛负 1
  opt <- try(nlminb(c(0.1, 0.1, 0.7), obj, lower=c(1e-6, 1e-6, 1e-6), upper=c(Inf, 1, 0.999)), silent=TRUE)
  if(inherits(opt, "try-error")) return(c(NA, NA, NA))
  return(opt$par)
}

est_GNGQMLE <- function(y, r = 1, is_symmetric = TRUE) {
  log_y2 <- 2 * log(abs(y) + 1e-150)
  theta_init <- if (r == 1) est_LQMLE(y) else est_GQMLE(y)
  if (any(is.na(theta_init))) return(c(NA, NA, NA)) 
  
  log_h_init <- garch_filter_log(theta_init[1], theta_init[2], theta_init[3], y)
  u_init <- log_y2 - log_h_init
  u_init <- pmin(u_init, 100)
  eta <- sign(y) * exp(0.5 * u_init) 
  
  best_ll <- -Inf
  best_family <- NULL
  best_shape <- NULL
  
  if (is_symmetric) {
    obj_gg <- function(p) {
      nu <- p[1]; lam <- p[2]
      if (nu <= 0 || lam <= 0) return(1e9)
      ll <- (nu / (2 * lam * gamma(1/nu))) * exp(-abs(eta/lam)^nu)
      -mean(log(ll + 1e-150)) 
    }
    opt_gg <- try(nlminb(c(1, 1), obj_gg, lower=c(0.3, 0.01), upper=c(10, 10)), silent=TRUE)
    if (!inherits(opt_gg, "try-error")) {
      best_ll <- -opt_gg$objective; best_family <- "gg"; best_shape <- opt_gg$par
    }
    
    obj_t <- function(p) {
      nu <- p[1]; lam <- p[2]
      if (nu <= 2 || lam <= 0) return(1e9)
      c_t <- gamma((nu+1)/2) / (lam * sqrt(nu*pi) * gamma(nu/2))
      ll <- log(c_t) - ((nu+1)/2) * log(1 + (eta^2)/(lam^2 * nu))
      -mean(ll)
    }
    opt_t <- try(nlminb(c(5, 1), obj_t, lower=c(2.01, 0.01), upper=c(50, 10)), silent=TRUE)
    if (!inherits(opt_t, "try-error")) {
      if (-opt_t$objective > best_ll) { 
        best_ll <- -opt_t$objective; best_family <- "t"; best_shape <- opt_t$par 
      }
    }
  } else {
    obj_sgg <- function(p) {
      beta <- p[1]; nu <- p[2]; lam <- p[3]
      if (beta <= 0 || abs(nu) >= 0.99 || lam <= 0) return(1e9)
      ll <- try(dsgg(eta, beta, nu, lam, log=TRUE), silent=TRUE)
      if (inherits(ll, "try-error") || any(is.na(ll))) return(1e9)
      -mean(ll)
    }
    opt_sgg <- try(nlminb(c(1, 0, 1), obj_sgg, lower=c(0.3, -0.99, 0.01), upper=c(5, 0.99, 10)), silent=TRUE)
    if (!inherits(opt_sgg, "try-error")) {
      best_ll <- -opt_sgg$objective; best_family <- "sgg"; best_shape <- opt_sgg$par
    }
    
    obj_piv <- function(p) {
      m <- p[1]; nu <- p[2]; lam <- p[3]
      if (m <= 0.5 || lam <= 0) return(1e9)
      ll <- try(dpearsonIV(eta, m=m, nu=nu, location=0, scale=lam, log=TRUE), silent=TRUE)
      if (inherits(ll, "try-error") || any(is.na(ll))) return(1e9)
      -mean(ll)
    }
    opt_piv <- try(nlminb(c(3, 0, 1), obj_piv, lower=c(0.51, -10, 0.01), upper=c(50, 10, 10)), silent=TRUE)
    if (!inherits(opt_piv, "try-error")) {
      if (-opt_piv$objective > best_ll) { 
        best_ll <- -opt_piv$objective; best_family <- "pIV"; best_shape <- opt_piv$par 
      }
    }
  }
  
  if (is.null(best_family)) return(theta_init) 
  
  final_obj <- function(theta) {
    w <- theta[1]; a <- theta[2]; b <- theta[3]
    if (w <= 0 || a < 0 || b < 0 || b >= 0.999) return(1e9)
    
    log_h <- garch_filter_log(w, a, b, y)
    u <- log_y2 - log_h
    u <- pmin(u, 100) 
    res <- sign(y) * exp(0.5 * u)
    
    if (best_family == "t") {
      nu <- best_shape[1]; lam <- best_shape[2]
      c_t <- gamma((nu+1)/2) / (lam * sqrt(nu*pi) * gamma(nu/2))
      ll <- log(c_t) - ((nu+1)/2) * log(1 + (res^2)/(lam^2 * nu))
    } else if (best_family == "gg") {
      nu <- best_shape[1]; lam <- best_shape[2]
      ll <- log((nu / (2 * lam * gamma(1/nu))) * exp(-abs(res/lam)^nu) + 1e-150)
    } else if (best_family == "sgg") {
      beta <- best_shape[1]; nu <- best_shape[2]; lam <- best_shape[3]
      ll <- dsgg(res, beta, nu, lam, log=TRUE)
    } else if (best_family == "pIV") {
      m <- best_shape[1]; nu <- best_shape[2]; lam <- best_shape[3]
      ll <- dpearsonIV(res, m=m, nu=nu, location=0, scale=lam, log=TRUE)
    }
    
    obj_val <- -0.5*log_h + ll
    -mean(obj_val)
  }
  
  # 宸蹭慨鏀? alpha 涓婄晫浠?5 鏀逛负 1
  opt_final <- try(nlminb(theta_init, final_obj, lower=c(1e-6, 1e-6, 1e-6), upper=c(Inf, 1, 0.999)), silent=TRUE)
  if (inherits(opt_final, "try-error")) return(theta_init)
  
  return(opt_final$par)
}

est_Fan2014 <- function(y) {
  log_y2 <- 2 * log(abs(y) + 1e-150)
  
  theta_init <- est_GQMLE(y)
  if (any(is.na(theta_init))) return(c(NA, NA, NA))
  
  log_h_init <- garch_filter_log(theta_init[1], theta_init[2], theta_init[3], y)
  u_init <- log_y2 - log_h_init
  eta_res <- sign(y) * soft_exp(0.5 * u_init, limit=150) 
  
  obj_eta <- function(eta) {
    if (eta <= 0) return(1e9)
    ll <- -4 * log(1 + (eta_res/eta)^2 / 7)
    -mean(-log(eta) + ll, na.rm=TRUE)
  }
  opt_eta <- try(nlminb(1, obj_eta, lower=0.01, upper=10), silent=TRUE)
  if (inherits(opt_eta, "try-error")) return(theta_init)
  eta_f <- opt_eta$par
  
  final_obj <- function(theta) {
    w <- theta[1]; a <- theta[2]; b <- theta[3]
    if (w <= 0 || a < 0 || b < 0 || b >= 0.999) return(1e9)
    
    log_h <- garch_filter_log(w, a, b, y)
    u <- log_y2 - log_h
    res <- sign(y) * soft_exp(0.5 * u, limit=150)
    
    ll <- -4 * log(1 + (res/eta_f)^2 / 7)
    -mean(-0.5 * log_h + ll, na.rm=TRUE)
  }
  
  # 宸蹭慨鏀? alpha 涓婄晫浠?5 鏀逛负 1
  opt_final <- try(nlminb(theta_init, final_obj, lower=c(1e-8, 1e-6, 1e-6), upper=c(Inf, 1, 0.999), scale=c(1, 1, 1)), silent=TRUE)
  if (inherits(opt_final, "try-error")) return(theta_init)
  
  return(opt_final$par)
}
# =========================================================================
# 5. 骞惰浠跨湡涓荤▼搴?
# =========================================================================
param_dict <- list(
  "stationary" = list(alpha=0.15, beta=0.6),
  "explosive"  = list(alpha=0.2,  beta=0.9),
  "boundary"   = list(
    "N(0,1)" = 0.109651, "t9" = 0.113516, "t7" = 0.115343, 
    "t5" = 0.120145, "t3" = 0.150827,
    "gg1" = 0.120694, "gg0.8" = 0.128205, "gg0.4" = 0.191652,
    "MixN" = 0.063321, "MixT" = 0.177295, 
    "pIV1" = 0.115693, "pIV2" = 0.120956, "sgg1" = 0.171513, "sgg2" = 0.206497
  )
)

num_cores <- detectCores() - 1
cl <- makeCluster(num_cores, outfile = "")
registerDoParallel(cl)
results_list <- list()

cat(sprintf("=== 浠跨湡寮€濮?| 鏍稿績鏁? %d | 鏃堕棿: %s ===\n\n", num_cores, Sys.time()))

for (dgp in dgps) {
  cat(sprintf("[%s] 姝ｅ湪璁＄畻 %s 鐨勭悊璁虹Н鍒嗗父鏁?..\n", Sys.time(), dgp))
  scales <- get_rescale_factors_exact(dgp)
  cat(sprintf("   => 绉垎瀹屾垚: C1=%.4f, C2=%.4f (鐞嗚涓?), C3=%.4f, RMS=%.4f\n", 
              scales$C1, scales$C2, scales$C3, scales$RMS))
  
  for (regime in regimes) {
    for (n in n_values) {
      cat(sprintf("[%s] 杩涘害: %s | %s | n=%d | 骞惰杩愮畻涓?..\n", 
                  Sys.time(), dgp, regime, n))
      
      w0 <- 0.25
      b0 <- if(regime == "stationary") 0.6 else 0.9
      a0 <- if(regime == "boundary") param_dict$boundary[[dgp]] else param_dict[[regime]]$alpha
      true_params <- c(w0, a0, b0)
      
      res_matrix <- foreach(rep = 1:replications, .combine = rbind, 
                            .export = c("rsgg", "dsgg", "rPE", "dgg", "rpearsonIV", "dpearsonIV", 
                                        "rmixN", "dmixN", "rmixT", "dmixT", 
                                        "generate_eta", "est_GNGQMLE", "est_Fan2014", "est_LQMLE", 
                                        "est_GQMLE", "est_LADE", "garch_filter_log", "soft_exp")) %dopar% {
                                          
                                          eta <- generate_eta(dgp, n + 500, scales$RMS) 
                                          h <- numeric(n + 500); y <- numeric(n + 500)
                                          h[1] <- w0 / (1 - a0 - b0)
                                          if(h[1] < 0 || is.na(h[1])) h[1] <- w0
                                          y[1] <- eta[1] * sqrt(h[1])
                                          
                                          for (t in 2:(n + 500)) {
                                            h[t] <- w0 + a0 * y[t-1]^2 + b0 * h[t-1]
                                            y[t] <- eta[t] * sqrt(h[t])
                                          }
                                          y_obs <- y[501:(n + 500)] 
                                          
                                          safe_est <- function(func) {
                                            res <- tryCatch(func(y_obs), error = function(e) rep(NA, 3))
                                            return(res)
                                          }
                                          
                                          # 灏嗚繖 4 涓柊澧炲垎甯冨姞鍏ュ绉伴樀钀ワ紒
                                          is_sym <- dgp %in% c("N(0,1)", "gg1", "t9", "t7", "t5", "t3", 
                                                               "gg0.8", "gg0.4", "MixN", "MixT")
                                          
                                          theta_gng1 <- safe_est(function(y) est_GNGQMLE(y, r=1, is_symmetric=is_sym))
                                          theta_gng2 <- safe_est(function(y) est_GNGQMLE(y, r=2, is_symmetric=is_sym))
                                          theta_fan  <- safe_est(function(y) est_Fan2014(y)) 
                                          theta_lq   <- safe_est(function(y) est_LQMLE(y))
                                          theta_gq   <- safe_est(function(y) est_GQMLE(y))
                                          theta_lade <- safe_est(function(y) est_LADE(y))
                                          
                                          rescale <- function(theta_star, c_scale) {
                                            if(any(is.na(theta_star))) return(rep(NA, 3))
                                            return(c(theta_star[1]/c_scale, theta_star[2]/c_scale, theta_star[3]))
                                          }
                                          
                                          c(rescale(theta_gng1, scales$C1), rescale(theta_gng2, scales$C2),
                                            rescale(theta_fan, scales$C2),
                                            rescale(theta_lq, scales$C1),   rescale(theta_gq, scales$C2),
                                            rescale(theta_lade, scales$C3))
                                        }
      
      calc_metrics <- function(est_matrix, true_val) {
        bias <- colMeans(est_matrix, na.rm = TRUE) - true_val
        rmse <- sqrt(colMeans((sweep(est_matrix, 2, true_val))^2, na.rm = TRUE))
        return(c(bias * 10, rmse * 10)) 
      }
      
      metrics_gng1 <- calc_metrics(res_matrix[, 1:3], true_params)
      metrics_gng2 <- calc_metrics(res_matrix[, 4:6], true_params)
      metrics_fan  <- calc_metrics(res_matrix[, 7:9], true_params)
      metrics_lq   <- calc_metrics(res_matrix[, 10:12], true_params)
      metrics_gq   <- calc_metrics(res_matrix[, 13:15], true_params)
      metrics_lade <- calc_metrics(res_matrix[, 16:18], true_params)
      
      row_df <- data.frame(
        DGP = dgp, Regime = regime, n = n, Metric = c("Bias", "RMSE"),
        GNG1_w=metrics_gng1[c(1,4)], GNG1_a=metrics_gng1[c(2,5)], GNG1_b=metrics_gng1[c(3,6)],
        GNG2_w=metrics_gng2[c(1,4)], GNG2_a=metrics_gng2[c(2,5)], GNG2_b=metrics_gng2[c(3,6)],
        Fan_w =metrics_fan[c(1,4)],  Fan_a =metrics_fan[c(2,5)],  Fan_b =metrics_fan[c(3,6)],
        LQ_w  =metrics_lq[c(1,4)],   LQ_a  =metrics_lq[c(2,5)],   LQ_b  =metrics_lq[c(3,6)],
        GQ_w  =metrics_gq[c(1,4)],   GQ_a  =metrics_gq[c(2,5)],   GQ_b  =metrics_gq[c(3,6)],
        LADE_w=metrics_lade[c(1,4)], LADE_a=metrics_lade[c(2,5)], LADE_b=metrics_lade[c(3,6)]
      )
      results_list[[length(results_list) + 1]] <- row_df
    }
  }
}
stopCluster(cl)

final_results <- do.call(rbind, results_list)
print(final_results)

split_results <- list(
  "Stationary" = subset(final_results, Regime == "stationary"),
  "Boundary"   = subset(final_results, Regime == "boundary"),
  "Explosive"  = subset(final_results, Regime == "explosive")
)

if (requireNamespace("writexl", quietly = TRUE)) {
  excel_name <- ifelse(TEST_MODE, "GNGQMLE_TEST_Results.xlsx", "GNGQMLE_Simulation_Results.xlsx")
  writexl::write_xlsx(split_results, path = excel_name)
  cat(sprintf("\n=== 浠跨湡缁撴潫 | 鏃堕棿: %s ===\n", Sys.time()))
  cat(sprintf("缁撴灉淇濆瓨鑷? %s\n", excel_name))
} else {
  prefix <- ifelse(TEST_MODE, "GNGQMLE_TEST", "GNGQMLE_Simulation")
  for (regime_name in names(split_results)) {
    sub_data <- split_results[[regime_name]]
    if(nrow(sub_data) > 0) {
      file_name <- sprintf("%s_%s.csv", prefix, regime_name)
      write.csv(sub_data, file_name, row.names = FALSE)
    }
  }
}