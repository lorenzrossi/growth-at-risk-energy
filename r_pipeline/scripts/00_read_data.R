# 00_read_data.R -- data construction.
# Replaces: ReadData.m, ReadData_OIL.m, Country_Read_Data.m
#
# Inputs (set paths in config.R):
#   - EU27 aggregate CSV with columns Time, IP, IP_t1, GAS_PRICE, GAS_PRICE_t1
#     (ue27_df.csv / DataIP_GaR.csv)
#   - Brent CSV with columns Time, BRENT_PRICE, BRENT_PRICE_t1 (DataOIL-style)
#   - Country/sector panel CSV with columns Time, nace_r2, geo, IP, IP_t1,
#     GAS_PRICE, GAS_PRICE_t1, BRENT_PRICE (22d0e721-style)
#   - Optional EU27 sector CSVs (e.g. MIG - energy) with Time, IP, GAS_PRICE
#
# Output: data/prepared.rds -- a named list of data.frames, each with
#   Time, Diff_IP, Diff_GP (or Diff_BP), NGPI / NOPI as relevant.

source(file.path(cfg$root, "R", "utils.R"))

read_time_csv <- function(path) {
  df <- utils::read.csv(path, check.names = FALSE)
  tcol <- intersect(c("Time", "Date", "TIME_PERIOD"), names(df))[1]
  df$Time <- as.Date(df[[tcol]])
  df
}

prep_series <- function(df, price_col = "GAS_PRICE", diff_name = "Diff_GP",
                        ni_name = "NGPI") {
  df <- df[order(df$Time), ]
  out <- data.frame(Time = df$Time)
  out$IP <- df$IP
  out$Diff_IP <- log_diff(df$IP)
  if (price_col %in% names(df)) {
    out[[price_col]] <- df[[price_col]]
    out[[diff_name]] <- log_diff(df[[price_col]])
    out[[ni_name]] <- net_increase(df[[price_col]])   # canonical 100*log construction
  }
  subsample_window(out)
}

data_list <- list()

## EU27 aggregate ------------------------------------------------------------
agg_raw <- read_time_csv(cfg$files$aggregate)
data_list$aggregate <- prep_series(agg_raw)

## Brent ---------------------------------------------------------------------
oil_raw <- read_time_csv(cfg$files$brent)
oil_raw <- oil_raw[order(oil_raw$Time), ]
oil <- data.frame(Time = oil_raw$Time,
                  BRENT_PRICE = oil_raw$BRENT_PRICE,
                  Diff_BP = log_diff(oil_raw$BRENT_PRICE),
                  NOPI = net_increase(oil_raw$BRENT_PRICE))
data_list$oil <- subsample_window(oil)

## Country x sector panel ----------------------------------------------------
panel <- read_time_csv(cfg$files$panel)
for (geo in unique(panel$geo)) {
  for (sec in unique(panel$nace_r2)) {
    sub <- panel[panel$geo == geo & panel$nace_r2 == sec, ]
    if (nrow(sub) < 24) next
    key <- paste(gsub("\\W+", "_", geo), gsub("\\W+", "_", sec), sep = ".")
    data_list[[key]] <- prep_series(sub)
  }
}

## Optional EU27 sector files (e.g. MIG energy) ------------------------------
for (nm in names(cfg$files$eu_sectors)) {
  path <- cfg$files$eu_sectors[[nm]]
  if (!file.exists(path)) {
    warning("EU sector file not found, skipping: ", path)
    next
  }
  sec_raw <- read_time_csv(path)
  data_list[[paste0("EU27.", nm)]] <- prep_series(sec_raw)
}

## Merge oil into every IP dataset by date (explicit join, not row order!) ---
for (nm in setdiff(names(data_list), "oil")) {
  d <- merge(data_list[[nm]],
             data_list$oil[, c("Time", "BRENT_PRICE", "Diff_BP", "NOPI")],
             by = "Time", all.x = TRUE)
  data_list[[nm]] <- d[order(d$Time), ]
}

dir.create(file.path(cfg$root, "data"), showWarnings = FALSE)
saveRDS(data_list, file.path(cfg$root, "data", "prepared.rds"))
message("00_read_data: prepared ", length(data_list), " datasets -> data/prepared.rds")
message("  keys: ", paste(names(data_list), collapse = ", "))
