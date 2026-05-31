## Create processed return series from the raw public data files.
## Run from the repository root with: Rscript R/create_processed_data.R

read_exchange_returns <- function(symbol, start_date, end_date) {
  path <- sprintf("data/raw/USD_%s_raw.csv", symbol)
  dat <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  names(dat)[1:2] <- c("date", "close")
  dat$date <- as.Date(dat$date, tryFormats = c("%Y-%m-%d", "%Y/%m/%d", "%Y-%m-%d"))
  dat$close <- as.numeric(dat$close)
  dat <- dat[is.finite(dat$close) & !is.na(dat$date), ]
  dat <- dat[order(dat$date), ]
  ret <- -diff(log(dat$close))
  date <- dat$date[-1]
  keep <- date >= as.Date(start_date) & date <= as.Date(end_date)
  data.frame(series = paste0("USD/", symbol), date = date[keep], return = as.numeric(ret[keep]))
}

read_crypto_returns <- function(symbol, start_date, end_date) {
  path <- sprintf("data/raw/%s-usd-max.csv", tolower(symbol))
  dat <- read.csv(path, stringsAsFactors = FALSE)
  dat$date <- as.Date(substr(dat$snapped_at, 1, 10))
  dat$price <- as.numeric(dat$price)
  dat <- dat[is.finite(dat$price) & !is.na(dat$date), ]
  dat <- dat[order(dat$date), ]
  ret <- diff(log(dat$price))
  date <- dat$date[-1]
  keep <- date >= as.Date(start_date) & date <= as.Date(end_date)
  data.frame(series = toupper(symbol), date = date[keep], return = as.numeric(ret[keep]))
}

exchange <- rbind(
  read_exchange_returns("TRY", "2021-01-01", "2024-12-31"),
  read_exchange_returns("ARS", "2016-01-01", "2019-12-31")
)
crypto <- rbind(
  read_crypto_returns("btc", "2017-01-01", "2020-12-31"),
  read_crypto_returns("eth", "2017-01-01", "2020-12-31"),
  read_crypto_returns("bnb", "2018-01-01", "2021-12-31"),
  read_crypto_returns("trx", "2018-01-01", "2021-12-31")
)

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
write.csv(exchange, "data/processed/exchange_returns.csv", row.names = FALSE)
write.csv(crypto, "data/processed/crypto_returns.csv", row.names = FALSE)
cat("Wrote processed exchange and cryptocurrency returns.\n")
