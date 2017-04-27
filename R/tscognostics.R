#' Compute a collection of cognostics for multivariate time series.
#' 
#' The function computes a set of cognostics to summarize statistical features
#' of multivariate time series.
#'
#' @param y a \code{ts} object. A multivariate time series.
#' @param normalise logical. If TRUE, scaling time series to be normally distributed
#' with mean 0 and 1.
#' @param width integer. A window size specified for variance change, level shift,
#' and lumpiness.
#'
#' @return Cognostics matrix for every time series.
#' 
#' @importFrom RcppRoll roll_mean
#' @importFrom RcppRoll roll_var
#' @importFrom ForeCA spectral_entropy
#' @importFrom mgcv gam
#' @export

tsmeasures <- function(y, normalise = TRUE, width) {
  y <- as.ts(y)
  tspy <- tsp(y)
  freq <- frequency(y)
  if (missing(width)) {
    width <- freq
  }
  if (width <= 1L) {
    stop("width should be more than 1.")
  }
  # Remove columns containing all NAs
  nay <- is.na(y)
  measures <- list()
  allna <- apply(nay, 2, all)
  x <- y[, !allna]
  measures$ACF1 <- apply(x, 2, FoAcf)
  if (normalise) {
    x <- as.ts(scale(x)) # Normalise data to mean = 0, sd = 1
  }
  trimx <- as.ts(apply(x, 2, Trim))
  tsp(trimx) <- tsp(x) <- tspy
  lumping <- apply(x, 2, Lump, width = width)
  measures$lumpiness <- sapply(lumping, function(x) x$lumpiness)
  measures$stationary <- sapply(lumping, function(x) x$stationary)
  measures$entropy <- apply(x, 2, Entropy)
  measures$lshift <- apply(trimx, 2, RLshift, width = width)
  measures$vchange <- apply(trimx, 2, RVarChange, width = width)
  measures$cpoints <- apply(x, 2, Cpoints)
  measures$fspots <- apply(x, 2, Fspots)
  missing <- apply(x, 2, CountNAs)
  measures$relativeNA <- sapply(missing, function(x) x$relativeNA)
  measures$fixedNA <- sapply(missing, function(x) x$fixedNA)
  # Mean and variance are not recommended for non-stationary time series
  # measures$mean <- colMeans(x, na.rm = TRUE)
  # measures$var <- apply(x, 2, var, na.rm = TRUE)
  varts <- apply(x, 2, VarTS, tspx = tspy)
  measures$trend <- sapply(varts, function(x) x$trend)
  measures$linearity <- sapply(varts, function(x) x$linearity)
  measures$curvature <- sapply(varts, function(x) x$curvature)
  measures$spikiness <- sapply(varts, function(x) x$spike)
  if (freq > 1) {
    measures$season <- sapply(varts, function(x) x$season)
    measures$peak <- sapply(varts, function(x) x$peak)
    measures$trough <- sapply(varts, function(x) x$trough)
  }
  tmp <- do.call(cbind, measures)
  nr <- ncol(y)
  nc <- length(measures)
  if (all(!allna)) {
    out <- tmp
  } else {
    mat <- matrix(, nrow = nr, ncol = nc)
    colnames(mat) <- colnames(tmp)
    mat[!allna, ] <- tmp
    out <- mat
  }
  return(out)
}

# Trimmed time series elimating outliers's influence
Trim <- function(x, trim = 0.1) {
  qtl <- quantile(x, c(trim, 1 - trim), na.rm = TRUE)
  lo <- qtl[1L]
  hi <- qtl[2L]
  x[x < lo | x > hi] <- NA
  return(x)
}

# Lumpiness and stationarity using sliding windows
Lump <- function(x, width) {
  nr <- length(x)
  lo <- seq(1, nr, by = width)
  up <- seq(width, nr + width, by = width)
  nsegs <- nr/width
  # Lumpiness
  varx <- sapply(1:nsegs, function(idx)
                 var(x[lo[idx]:up[idx]], na.rm = TRUE))
  lumpiness <- var(varx, na.rm = TRUE)
  # Degree of stationary
  meanx <- sapply(1:nsegs, function(idx) mean(x[lo[idx]:up[idx]], na.rm = TRUE))
  stationary <- var(meanx, na.rm = TRUE)
  return(list(lumpiness = lumpiness, stationary = stationary))
}

# Spectral entroy (using spectral_entropy from ForeCA package)
Entropy <- function(x) {
  entropy <- try(spectral_entropy(na.contiguous(x))[1L], silent = TRUE)
  if (class(entropy) == "try-error") {
    entropy <- NA
  }
  return(entropy)
}

# Autocorrelation function at lag 1
FoAcf <- function(x) {
  return(acf(x, plot = FALSE, na.action = na.exclude)$acf[2L])
}

# Level shift using rolling window
RLshift <- function(x, width) {
  rollmean <- roll_mean(x, width, na.rm = TRUE)
  lshifts <- tryCatch(max(abs(diff(rollmean, width)), na.rm = TRUE),
                      warning = function(w) w)
  if (any(class(lshifts) == "warning")) {
    lshifts <- NA
  }
  return(lshifts)
}

RVarChange <- function(x, width) {
  rollvar <- roll_var(x, width, na.rm = TRUE)
  vchange <- tryCatch(max(abs(diff(rollvar, width)), na.rm = TRUE),
                      warning = function(w) w)
  if (any(class(vchange) == "warning")) {
    vchange <- NA
  }
  return(vchange)
}

# Number of crossing points
Cpoints <- function(x) {
  nax <- is.na(x)
  starting <- which.min(nax)
  midline <- diff(range(x, na.rm = TRUE))/2
  ab <- x <= midline
  lenx <- length(x)
  p1 <- ab[1:(lenx - 1)]
  p2 <- ab[2:lenx]
  cross <- (p1 & !p2) | (p2 & !p1)
  out <- sum(cross, na.rm = TRUE)/(lenx - starting + 1)
  return(out)
}

# Flat spots using discretization
Fspots <- function(x) {
  cutx <- try(cut(x, breaks = 10, include.lowest = TRUE, labels = FALSE),
              silent = TRUE)
  if (class(cutx) == "try-error") {
    fspots <- NA
  } else {
    rlex <- rle(cutx)
    # Any flat spot
    return(max(rlex$lengths))
    # Low flat spots
    # ones <- (rlex$values == 1)
    # return(max(rlex$lengths[ones]))
  }
}

# Strength of trend and seasonality and spike
VarTS <- function(x, tspx) {
  x <- as.ts(x)
  tsp(x) <- tspx
  freq <- tspx[3]
  contx <- try(na.contiguous(x), silent = TRUE)
  len.contx <- length(contx)
  if (length(contx) <= 2 * freq || class(contx) == "try-error") {
    trend <- linearity <- curvature <- season <- spike <- peak <- trough <- NA
  } else {
    if (freq > 1L) {
      all.stl <- stl(contx, s.window = "periodic", robust = TRUE)
      starty <- start(contx)[2L]
      pk <- (starty + which.max(all.stl$time.series[, "seasonal"]) - 1L) %% freq
      th <- (starty + which.min(all.stl$time.series[, "seasonal"]) - 1L) %% freq
      pk <- ifelse(pk == 0, freq, pk)
      th <- ifelse(th == 0, freq, th)
      trend0 <- all.stl$time.series[, "trend"]
      fits <- trend0 + all.stl$time.series[, "seasonal"]
      adj.x <- contx - fits
      v.adj <- var(adj.x, na.rm = TRUE)
      detrend <- contx - trend0
      deseason <- contx - all.stl$time.series[, "seasonal"]
      peak <- pk * max(all.stl$time.series[, "seasonal"], na.rm = TRUE)
      trough <- th * min(all.stl$time.series[, "seasonal"], na.rm = TRUE)
      remainder <- all.stl$time.series[, "remainder"]
      season <- ifelse(var(detrend, na.rm = TRUE) < 1e-10, 0,
                       max(0, min(1, 1 - v.adj/var(detrend, na.rm = TRUE))))
    } else { # No seasonal component
      tt <- 1:len.contx
      trend0 <- fitted(mgcv::gam(contx ~ s(tt)))
      remainder <- contx - trend0
      deseason <- contx - trend0
      v.adj <- var(trend0, na.rm = TRUE)
    }
    trend <- ifelse(var(deseason, na.rm = TRUE) < 1e-10, 0,
                    max(0, min(1, 1 - v.adj/var(deseason, na.rm = TRUE))))
    n <- length(remainder)
    v <- var(remainder, na.rm = TRUE)
    d <- (remainder - mean(remainder, na.rm = TRUE))^2
    varloo <- (v * (n - 1) - d)/(n - 2)
    spike <- var(varloo, na.rm = TRUE)
    pl <- poly(1:len.contx, degree = 2)
    tren.coef <- coef(lm(trend0 ~ pl))[2:3]
    linearity <- tren.coef[1]
    curvature <- tren.coef[2]
  }
  if (freq > 1) { 
    return(list(trend = trend, season = season, spike = spike,
                peak = peak, trough = trough, linearity = linearity,
                curvature = curvature))
  } else { # No seasonal component
    return(list(trend = trend, spike = spike, linearity = linearity,
                curvature = curvature))
  }
}

# Count missing values
CountNAs <- function(x) {
  nonnax <- is.na(x)
  idx <- which(nonnax)
  starting <- idx[1L]
  ending <- idx[length(idx)]
  count_nas <- sum(nonnax)
  out1 <- count_nas/(ending - starting + 1)
  out2 <- count_nas/length(x)
  return(list(relativeNA = out1, fixedNA = out2))
}
