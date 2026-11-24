#!/usr/bin/env Rscript
## =====================================================================
## RGrowthCurves  v0.6.0
## Yeast growth rate (mu), doubling time, lag and area under the curve
## from multi-well plate-reader tables.
##
## Input : "vertical" table, one column per well, one row per read.
##         Header row = well/condition names. Optional first column
##         with time (header "Time", "t", "Vremya" (Cyrillic); minutes, hours or h:mm:ss).
##         Separator ; , or TAB and decimal point or comma are detected
##         automatically. UTF-8 (with or without BOM) or Windows-1251.
## Output: <input>_report.csv  - one row per well (AUC in OD*h, trapezoid rule)
##         plots/*.png         - one diagnostic plot per well
##         <input>_overview_NN.png - 4x4 overview pages with the fits
##
## Usage : Rscript RGrowthCurves.R input.csv [interval_min]
##         or edit SETTINGS and source() the file in RStudio.
## Base R only, no packages required.
## =====================================================================

VERSION <- "0.6.0"
LN2     <- log(2)

## ----------------------------- SETTINGS ------------------------------
## All user-tunable parameters live in this one list; main() passes it
## down to every function, so nothing depends on hidden global state.
CONFIG <- list(
  input_file        = "inn.csv",
  interval_min      = 5,          # minutes between reads (ignored if file has a time column)
  window_min        = 150,        # sliding-window width for the main Mu, min
  extra_windows_min = c(250),     # Mu is also reported for these windows (v0.4 used 250 min),
                                  # columns "Mu_w250" etc. - for comparison with older runs
  first_fit_min     = 0,          # windows may not start before this time, min
  min_r2            = 0.95,       # windows with lower R^2 are ignored when searching max slope
  blank_od          = 0,          # known background OD subtracted before log2 and AUC (0 = off).
                                  # Leave 0 if the plate has no empty-medium wells.
  baseline_min      = 120,        # net AUC baseline (plate + medium + inoculum), no blanks needed:
                                  # lowest smoothed OD from 0 to the start of the exponential-phase
                                  # fit (at least this many min; used alone if no fit)
  start_od_range    = c(5, 35),   # time range (min) averaged as starting OD (v0.4: rows 2-8)
  od_at_h           = 10,         # time point (h) for the OD10h columns
  od_at_half_min    = 30,         # +/- minutes averaged around od_at_h
  mu_cutoffs_h      = c(10, 12, 14),  # Mu_10 etc: max slope restricted by time, see next line
  mu_cutoff_rule    = "end",      # "end":   only windows that END before the cutoff, i.e. Mu as if
                                  #          the plate had been stopped at 10/12/14 h (comparable)
                                  # "start": windows STARTING before cutoff (v0.4 behaviour; such
                                  #          windows use data up to cutoff + window width)
  auc_cutoff_h      = 10,         # AUC is also reported for 0..auc_cutoff_h
  min_delta_od      = 0.05,       # wells growing less than this are flagged "no growth"
  out_dir           = "RGrowthCurves_output",
  out_sep           = ";",        # report separator
  out_dec           = ".",        # report decimal mark ("," for Russian-locale Excel)
  make_plots        = TRUE,
  overview_grid     = c(4, 4)     # rows x columns of wells per overview page
)
## ---------------------------------------------------------------------

## ============================ INPUT ==================================

## Read a text file that may be UTF-8 (with/without BOM) or Windows-1251
## (Excel on Russian Windows); drop blank lines.
read_lines_any_encoding <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  if (!all(validUTF8(lines))) {
    lines <- iconv(lines, from = "CP1251", to = "UTF-8")
  }
  lines[1] <- sub("^\ufeff", "", lines[1])          # strip BOM
  lines[nzchar(trimws(lines))]
}

## Guess the column separator: ";", TAB, else "," (a comma is trusted
## only when the other two are absent, since it may be a decimal mark).
detect_separator <- function(lines) {
  probe <- head(lines, 30)
  n_fields <- function(sep) median(lengths(strsplit(probe, sep, fixed = TRUE)))
  counts <- c(";" = n_fields(";"), "\t" = n_fields("\t"))
  if (max(counts) >= 2) return(names(counts)[which.max(counts)])
  ","
}

## Lenient string -> number: strips (non-breaking) spaces, decimal comma -> point.
to_number <- function(x) {
  x <- gsub("[[:space:]\u00a0]", "", x)
  x <- sub(",", ".", x, fixed = TRUE)
  suppressWarnings(as.numeric(x))
}

## Parse a time column to minutes. Accepts h:mm[:ss] or plain numbers;
## plain numbers are scaled by units mentioned in the header
## (hours/chas -> *60, seconds -> /60, minutes as is).
parse_time_column <- function(x, header) {
  x <- trimws(x)
  if (any(grepl(":", x, fixed = TRUE))) {           # h:mm:ss or h:mm
    parts <- lapply(strsplit(x, ":", fixed = TRUE), to_number)
    return(vapply(parts, function(v) {
      v <- c(v, rep(0, 3 - length(v)))[1:3]
      v[1] * 60 + v[2] + v[3] / 60
    }, numeric(1)))
  }
  v <- to_number(x)
  if (grepl("(?i)(hour|\\(h\\)|, *h\\b|\u0447\u0430\u0441|\\(\u0447\\))", header, perl = TRUE)) v <- v * 60
  if (grepl("(?i)(sec|\\(s\\)|, *s\\b|\u0441\u0435\u043a)", header, perl = TRUE))              v <- v / 60
  v
}

## Read the plate-reader table. Returns list(time = minutes per read,
## od = data frame with one numeric column per well).
read_growth_table <- function(path, interval_min) {
  lines <- read_lines_any_encoding(path)
  sep   <- detect_separator(lines)
  raw <- read.table(text = lines, sep = sep, header = TRUE, quote = "\"",
                    colClasses = "character", check.names = FALSE,
                    comment.char = "", strip.white = TRUE,
                    na.strings = c("", "NA", "OVRFLW", "OVER", "-", "*"),
                    encoding = "UTF-8", fill = TRUE)

  header   <- names(raw)
  has_time <- grepl("(?i)^\\s*(time|t|\u0432\u0440\u0435\u043c\u044f|\u0447\u0430\u0441)(\\s|,|\\(|_|$)",
                    header[1], perl = TRUE)
  if (has_time) {
    time_min <- parse_time_column(raw[[1]], header[1])
    raw      <- raw[-1]
  } else {
    time_min <- interval_min * (seq_len(nrow(raw)) - 1)
  }

  od <- as.data.frame(lapply(raw, to_number), check.names = FALSE)
  names(od) <- trimws(names(raw))
  ## drop empty columns (trailing separators) and rows without any data
  keep_col <- colSums(!is.na(od)) > 0 & nzchar(trimws(names(od)))
  keep_row <- rowSums(!is.na(od)) > 0 & !is.na(time_min)
  od       <- od[keep_row, keep_col, drop = FALSE]
  time_min <- time_min[keep_row]

  message(sprintf("Read %s: separator '%s', %d wells x %d reads, %s, %.1f h",
                  basename(path), ifelse(sep == "\t", "TAB", sep), ncol(od), nrow(od),
                  if (has_time) "time column found" else
                    sprintf("time from interval %g min", interval_min),
                  max(time_min) / 60))
  list(time = time_min, od = od)
}

## ============================ ANALYSIS ===============================

## Least-squares line through (t, y) without lm(): ~100x faster, same result.
## Returns c(slope, intercept, r2); NAs if fewer than 3 finite points.
fit_line <- function(t, y) {
  ok <- is.finite(y)
  if (sum(ok) < 3) return(c(slope = NA_real_, intercept = NA_real_, r2 = NA_real_))
  t <- t[ok]; y <- y[ok]
  var_t <- var(t); var_y <- var(y)
  slope <- cov(t, y) / var_t
  c(slope     = slope,
    intercept = mean(y) - slope * mean(t),
    r2        = if (var_y > 0) cov(t, y)^2 / (var_t * var_y) else NA_real_)
}

## Trapezoidal area under y(t) between t_from and t_to (min), in y*h.
## Endpoints are linearly interpolated; NA if the run does not reach t_to.
auc_trapezoid <- function(t, y, t_from = min(t), t_to = max(t)) {
  ok <- is.finite(y); t <- t[ok]; y <- y[ok]
  if (length(t) < 2 || t_to > max(t) + 1e-9) return(NA_real_)
  edges  <- approx(t, y, xout = c(t_from, t_to))$y   # interpolate at the limits
  inside <- t > t_from & t < t_to
  tt <- c(t_from, t[inside], t_to)
  yy <- c(edges[1], y[inside], edges[2])
  sum(diff(tt) * (head(yy, -1) + tail(yy, -1)) / 2) / 60
}

## Mean of y over t in [from, to] (NA if all missing).
mean_between <- function(t, y, from, to) {
  v <- y[t >= from & t <= to]
  if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)
}

## Fit a line to log2 OD in every sliding window of width window_min.
## Returns a data frame with one row per window: start, end, slope, intercept, r2.
fit_sliding_windows <- function(t, log2_od, window_min, first_fit_min) {
  ## window width in points: w intervals = w+1 reads
  w      <- max(2, round(window_min / median(diff(t))))
  starts <- which(t >= first_fit_min & seq_along(t) + w <= length(t))
  if (!length(starts))
    return(data.frame(start = numeric(0), end = numeric(0), slope = numeric(0),
                      intercept = numeric(0), r2 = numeric(0)))
  fits <- t(vapply(starts, function(i) fit_line(t[i:(i + w)], log2_od[i:(i + w)]),
                   c(slope = 0, intercept = 0, r2 = 0)))
  data.frame(start = t[starts], end = t[starts + w], fits)
}

## A window counts as exponential growth only if its slope is positive
## and the linear fit is good enough (R^2 >= min_r2).
is_usable_window <- function(windows, min_r2) {
  with(windows, is.finite(slope) & slope > 0 & !is.na(r2) & r2 >= min_r2)
}

## A well "grew" if its high quantile rises at least min_delta_od above
## the starting OD (quantile is robust to single spikes).
well_grew <- function(od, start_od, min_delta_od) {
  rise <- quantile(od, 0.95, na.rm = TRUE) - start_od
  is.finite(rise) && rise >= min_delta_od
}

## The usable window with the steepest slope = the exponential phase.
## Returns index plus the window's fit; all-NA fields when nothing is usable.
pick_best_window <- function(windows, usable) {
  empty <- list(index = NA_integer_, slope = NA_real_, intercept = NA_real_,
                r2 = NA_real_, start = NA_real_, end = NA_real_)
  if (!any(usable)) return(empty)
  idx <- which(usable)[which.max(windows$slope[usable])]
  c(list(index = idx), as.list(windows[idx, ]))
}

## Human-readable QC notes; "ok" when there is nothing to report.
quality_flags <- function(t, od_raw, usable, best, growth, cfg) {
  flags <- character(0)
  if (!growth) {
    flags <- c(flags, sprintf("no growth (dOD<%g)", cfg$min_delta_od))
  } else if (!any(usable)) {
    flags <- c(flags, sprintf("no window with R2>=%g", cfg$min_r2))
  }
  if (!is.na(best$index) &&
      max(od_raw[t >= best$start & t <= best$end], na.rm = TRUE) > 1.0)
    flags <- c(flags, "fit window reaches OD>1, check reader linearity")
  dt <- median(diff(t))
  if (max(t) < max(cfg$mu_cutoffs_h) * 60 - dt / 2)
    flags <- c(flags, sprintf("run shorter than %g h", max(cfg$mu_cutoffs_h)))
  if (length(flags)) paste(flags, collapse = "; ") else "ok"
}

## Max slope (as Mu, 1/h) among windows that lie entirely before each
## cutoff ("end" rule) or start before it ("start" rule) - i.e. the Mu
## one would get if the plate had been stopped at 10/12/14 h.
## Unlike the main Mu, these columns skip the R^2 and growth filters so a
## value is always reported; the matching R2_* columns show the fit
## quality of the selected window so the reader can judge it. NA only
## when the run is too short for the cutoff or no positive-slope window
## fits before it (e.g. a flat or declining curve).
mu_before_cutoffs <- function(t, windows, cfg) {
  dt       <- median(diff(t))
  positive <- is.finite(windows$slope) & windows$slope > 0
  fits <- vapply(cfg$mu_cutoffs_h, function(h) {
    if (max(t) < h * 60 - dt / 2) return(c(mu = NA_real_, r2 = NA_real_))  # run too short
    edge <- if (cfg$mu_cutoff_rule == "end") windows$end else windows$start
    sel  <- positive & edge <= h * 60
    if (!any(sel)) return(c(mu = NA_real_, r2 = NA_real_))
    idx <- which(sel)[which.max(windows$slope[sel])]
    c(mu = LN2 * 60 * windows$slope[idx], r2 = windows$r2[idx])
  }, c(mu = 0, r2 = 0))
  list(mu = fits["mu", ], r2 = fits["r2", ])
}

## Mu (1/h) for additional window widths, for comparison with older runs.
mu_with_extra_windows <- function(t, log2_od, growth, cfg) {
  vapply(cfg$extra_windows_min, function(window_min) {
    windows <- fit_sliding_windows(t, log2_od, window_min, cfg$first_fit_min)
    usable  <- is_usable_window(windows, cfg$min_r2) & growth
    if (any(usable)) LN2 * 60 * max(windows$slope[usable]) else NA_real_
  }, numeric(1))
}

## Baseline for net AUC: lowest OD early in the run, after a running
## median of 5 reads which removes single-read spikes (bubbles,
## condensation). "Early" = up to the start of the exponential-phase
## fit, but always at least baseline_min minutes.
early_baseline <- function(t, od, fit_start, baseline_min) {
  ok <- !is.na(od)
  smoothed <- od
  if (sum(ok) >= 5) smoothed[ok] <- runmed(od[ok], 5, endrule = "median")
  early <- ok & t <= max(baseline_min, fit_start, na.rm = TRUE)
  if (any(early)) min(smoothed[early]) else NA_real_
}

## Full analysis of one well. Returns a list of report values (mu, tau,
## lag, AUCs, qc, ...) plus the internals needed by the diagnostic plots
## (windows, usable, best, log2_od, od).
analyse_well <- function(t, od_raw, cfg) {
  od      <- od_raw - cfg$blank_od
  log2_od <- suppressWarnings(log2(od))   # OD <= 0 -> NaN/-Inf, skipped in fits

  start_od <- mean_between(t, od, cfg$start_od_range[1], cfg$start_od_range[2])
  growth   <- well_grew(od, start_od, cfg$min_delta_od)

  windows <- fit_sliding_windows(t, log2_od, cfg$window_min, cfg$first_fit_min)
  usable  <- is_usable_window(windows, cfg$min_r2) & growth
  best    <- pick_best_window(windows, usable)

  qc       <- quality_flags(t, od_raw, usable, best, growth, cfg)
  mu_cut   <- mu_before_cutoffs(t, windows, cfg)
  mu_extra <- mu_with_extra_windows(t, log2_od, growth, cfg)
  baseline <- early_baseline(t, od, best$start, cfg$baseline_min)

  t_cut <- cfg$auc_cutoff_h * 60
  od_at <- mean_between(t, od,
                        cfg$od_at_h * 60 - cfg$od_at_half_min,
                        cfg$od_at_h * 60 + cfg$od_at_half_min)

  list(
    mu          = LN2 * 60 * best$slope,                    # specific growth rate, 1/h
    tau         = 1 / (60 * best$slope),                    # doubling time, h
    start_od    = start_od,
    od_at_minus = od_at - start_od,
    od_at       = od_at,
    win_start   = best$start,                               # start of max-slope window, min
    lag_corr    = (log2(start_od) - best$intercept) / best$slope,  # tangent lag, min
    r2          = best$r2,
    mu_cut      = mu_cut$mu,
    r2_cut      = mu_cut$r2,
    mu_extra    = mu_extra,
    baseline    = baseline,
    auc         = auc_trapezoid(t, od),                     # OD*h, whole run
    auc_cut     = auc_trapezoid(t, od, min(t), t_cut),      # OD*h, 0..cutoff (comparable)
    auc_net     = auc_trapezoid(t, od - baseline),
    auc_cut_net = auc_trapezoid(t, od - baseline, min(t), t_cut),
    qc          = qc,
    ## internals for the plots
    windows = windows, usable = usable, best = best,
    log2_od = log2_od, od = od
  )
}

## ============================ PLOTS ==================================

## Filesystem-safe well name for the per-well PNG files.
make_safe_filename <- function(x) {
  x <- iconv(x, "UTF-8", "ASCII//TRANSLIT", sub = "")
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  substr(gsub("^_|_$", "", x), 1, 40)
}

## One diagnostic plot per well: raw OD, log2 OD, every window's slope
## and R^2, the best exponential-phase fit and the tangent lag.
plot_well <- function(file, title, t, res) {
  png(file, width = 1000, height = 700, res = 110)
  on.exit(dev.off())
  finite_log <- res$log2_od[is.finite(res$log2_od)]
  y_min <- if (length(finite_log)) min(-4, floor(min(finite_log))) else -4
  has_fit <- !is.na(res$best$index)

  plot(t, res$log2_od, ylim = c(y_min, 2), col = "blue", pch = 1, cex = 0.6,
       main = title, xlab = "Time, min", ylab = "OD (black) / log2 OD (blue)")
  if (has_fit)
    rect(res$best$start, y_min - 1, res$best$end, 3,
         col = adjustcolor("red", 0.07), border = NA)
  points(t, res$od, pch = 1, cex = 0.6)

  ok <- is.finite(res$windows$slope)
  points(res$windows$start[ok], 10 * LN2 * 60 * res$windows$slope[ok] + y_min,  # Mu x10, offset to bottom
         col = ifelse(res$usable[ok], "darkgreen", "grey70"), pch = 20, cex = 0.5)
  points(res$windows$start[ok], res$windows$r2[ok], col = "orange", pch = 20, cex = 0.5)

  if (is.finite(log2(res$start_od))) abline(h = log2(res$start_od), lty = 2, col = "darkgreen")
  if (has_fit) {
    abline(a = res$best$intercept, b = res$best$slope, col = "red")
    abline(v = res$win_start, col = "blue")
    points(res$lag_corr, log2(res$start_od), col = "red", pch = 19, cex = 1.5)
    legend("topleft", bty = "n", text.col = "red", cex = 1.3,
           legend = c(sprintf("mu = %.3f 1/h", res$mu),
                      sprintf("tau = %.2f h", res$tau),
                      sprintf("lag = %.0f min", res$lag_corr)))
  } else {
    legend("topleft", bty = "n", text.col = "red", legend = res$qc)
  }
  legend("bottomright", bty = "n", cex = 0.8,
         pch = c(1, 1, 20, 20, NA, NA), lty = c(NA, NA, NA, NA, 1, 2),
         col = c("black", "blue", "darkgreen", "orange", "red", "darkgreen"),
         legend = c("OD", "log2 OD", sprintf("10*mu%+d (window start)", y_min),
                    "R2 of window", "best fit", "log2 starting OD"))
}

## 4x4 overview pages: raw OD curves with the fitted exponential segment,
## all panels on the same y scale so wells are visually comparable.
plot_overview <- function(stem, t, od, results, cfg) {
  wells_per_page <- prod(cfg$overview_grid)
  pages <- split(seq_len(ncol(od)), ceiling(seq_len(ncol(od)) / wells_per_page))
  y_max <- max(od, na.rm = TRUE)
  for (p in seq_along(pages)) {
    png(file.path(cfg$out_dir, sprintf("%s_overview_%02d.png", stem, p)),
        width = 420 * cfg$overview_grid[2], height = 330 * cfg$overview_grid[1], res = 100)
    par(mfrow = cfg$overview_grid, mar = c(3, 3, 2.2, 0.6), mgp = c(1.7, 0.5, 0),
        oma = c(0, 0, 2, 0))
    for (j in pages[[p]]) {
      res <- results[[j]]
      plot(t / 60, od[[j]], type = "n", ylim = c(0, y_max), xlab = "Time, h", ylab = "OD",
           main = sprintf("%d. %s", j, substr(names(od)[j], 1, 42)),
           cex.main = 0.9, font.main = 1)
      if (!is.na(res$best$index)) {
        rect(res$best$start / 60, -1, res$best$end / 60, 2 * y_max,
             col = adjustcolor("red", 0.08), border = NA)
        tt <- seq(res$best$start, res$best$end, length.out = 40)
        lines(tt / 60, 2^(res$best$intercept + res$best$slope * tt) + cfg$blank_od,
              col = "red", lwd = 3)
      }
      lines(t / 60, od[[j]], col = if (res$qc == "ok") "black" else "grey55")
      if (!is.na(res$best$index))
        legend("topleft", bty = "n", text.col = "red", cex = 0.95,
               legend = sprintf("mu = %.3f  tau = %.2f h", res$mu, res$tau))
      if (res$qc != "ok")
        legend("bottomright", bty = "n", text.col = "grey30", cex = 0.75, legend = res$qc)
    }
    mtext(sprintf("%s   page %d/%d   (window %g min, red = max-slope fit)",
                  stem, p, length(pages), cfg$window_min), outer = TRUE, line = 0.5)
    dev.off()
  }
}

## ============================ REPORT =================================

## Assemble the per-well results into the report table.
## Column names are kept compatible with v0.4; new columns go at the end.
build_report <- function(results, well_names, cfg) {
  num_label     <- function(x) format(x, trim = TRUE)   # "10" not "10.0" in column names
  result_column <- function(name) sapply(results, `[[`, name)

  report <- data.frame(conditions = well_names,
                       Mu = result_column("mu"),
                       check.names = FALSE)
  report[["tau12, h"]]        <- result_column("tau")
  report[["starting OD"]]     <- result_column("start_od")
  report[[sprintf("OD%sh_minusStartingOD", num_label(cfg$od_at_h))]] <- result_column("od_at_minus")
  report[[sprintf("OD%sh_abs", num_label(cfg$od_at_h))]]             <- result_column("od_at")
  report[["lag, min"]]            <- result_column("win_start")
  report[["lag_corrected, min"]]  <- result_column("lag_corr")
  report[["r2"]]                  <- result_column("r2")
  report[["RGScript_version"]]    <- VERSION

  ## one Mu + R2 column pair per time cutoff, one Mu column per extra window width
  mu_cut <- do.call(rbind, lapply(results, `[[`, "mu_cut"))
  r2_cut <- do.call(rbind, lapply(results, `[[`, "r2_cut"))
  for (k in seq_along(cfg$mu_cutoffs_h)) {
    label <- num_label(cfg$mu_cutoffs_h[k])
    report[[paste0("Mu_", label)]] <- mu_cut[, k]
    report[[paste0("R2_", label)]] <- r2_cut[, k]
  }
  mu_extra <- do.call(rbind, lapply(results, `[[`, "mu_extra"))
  colnames(mu_extra) <- paste0("Mu_w", num_label(cfg$extra_windows_min))
  report <- cbind(report, mu_extra)

  report[["AUC, OD*h"]]        <- result_column("auc")
  report[[sprintf("AUC_0-%sh, OD*h", num_label(cfg$auc_cutoff_h))]]     <- result_column("auc_cut")
  report[["AUC_net, OD*h"]]    <- result_column("auc_net")
  report[[sprintf("AUC_0-%sh_net, OD*h", num_label(cfg$auc_cutoff_h))]] <- result_column("auc_cut_net")
  report[["baseline OD"]]      <- result_column("baseline")
  report[["fit window, min"]]  <- cfg$window_min
  report[["Mu_cutoff_rule"]]   <- cfg$mu_cutoff_rule
  report[["QC"]]               <- result_column("qc")

  is_numeric <- vapply(report, is.numeric, logical(1))
  report[is_numeric] <- lapply(report[is_numeric], signif, digits = 5)
  report
}

## Write the report with a BOM so Excel shows Cyrillic/Greek names correctly.
write_report <- function(report, path, cfg) {
  con <- file(path, open = "w", encoding = "UTF-8")
  cat("\ufeff", file = con)
  write.table(report, con, sep = cfg$out_sep, dec = cfg$out_dec,
              row.names = FALSE, na = "")
  close(con)
}

## ============================ MAIN ===================================

## Command line overrides: Rscript RGrowthCurves.R input.csv [interval_min]
apply_cli_args <- function(cfg) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1) cfg$input_file   <- args[1]
  if (length(args) >= 2) cfg$interval_min <- as.numeric(sub(",", ".", args[2], fixed = TRUE))
  cfg
}

main <- function() {
  cfg  <- apply_cli_args(CONFIG)
  data <- read_growth_table(cfg$input_file, cfg$interval_min)

  dir.create(cfg$out_dir, showWarnings = FALSE)
  if (cfg$make_plots) dir.create(file.path(cfg$out_dir, "plots"), showWarnings = FALSE)
  stem <- sub("\\.[^.]*$", "", basename(cfg$input_file))

  results <- vector("list", ncol(data$od))
  for (j in seq_along(data$od)) {
    title <- names(data$od)[j]
    res   <- analyse_well(data$time, data$od[[j]], cfg)
    results[[j]] <- res
    if (cfg$make_plots)
      plot_well(file.path(cfg$out_dir, "plots",
                          sprintf("%03d_%s.png", j, make_safe_filename(title))),
                title, data$time, res)
    message(sprintf("%3d  %-45s mu=%6.3f  tau=%5.2f h  %s",
                    j, substr(title, 1, 45), res$mu, res$tau, res$qc))
  }

  report   <- build_report(results, names(data$od), cfg)
  out_file <- file.path(cfg$out_dir, paste0(stem, "_report.csv"))
  write_report(report, out_file, cfg)
  if (cfg$make_plots) plot_overview(stem, data$time, data$od, results, cfg)
  message("Report written to ", out_file)
}

main()
