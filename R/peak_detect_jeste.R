peak_detect_jeste <- function(x, y, span.range = c(0.25, 1),
                              span.inc = 0.05, k.fold = 5,
                              find.peak = TRUE, find.dip = TRUE,
                              peak.range = c(0.5, 0.99),
                              dip.range = c(0.01, 0.5)) {
    # Peak or dip ranges must not lie on the boundary of the ERP time domain
    if ((peak.range[2] >= max(x)) |
        (peak.range[1] <= min(x)) |
        (dip.range[2] >= max(x)) |
        (dip.range[1] <= min(x))) {
        stop("Feature ranges must not contain endpoints of ERP time domain")
    }
    
    # Take cross-sectional mean of subset of X_ijkl(t)
    # mean.time <- apply(sub.erp, 2, mean)
    # Implement loess smoothing with cross-validation
    loess_smooth <- loess.wrapper(x, y,
                                  span.vals=seq(span.range[1],
                                                span.range[2], span.inc),
                                  folds=k.fold)
    loess.smooth <- list(x=loess_smooth$x, y=loess_smooth$fitted)
    # Peak feature
    if (find.peak) {
        # Indicator for boundary issues. If peak estimate is located on
        # twice the size of the original feature boundary, the estimate
        # will be considered missing
        boundary <- 0
        # Interval for identifying peak
        peak.interval <- which(loess.smooth$x >= peak.range[1] &
                                   loess.smooth$x <= peak.range[2])
        # Check if peak is on the boundary and iterate
        repeat {
            # Identify maximum within interval
            peak.ind <- which.max(loess.smooth$y[peak.interval])
            # Shift interval if maximum on boundary
            if (peak.ind == length(peak.interval)) {
                peak.interval <- peak.interval + 1
            } else if (peak.ind == 1) {
                peak.interval <- peak.interval - 1
            } else {
                break
            } # end boundary if statement for peak
            # Quality check to see if estimate is far from original interval
            half.range <- (peak.range[2] - peak.range[1])/2
            if ((max(x[peak.interval]) > peak.range[2] + half.range) |
                (min(x[peak.interval]) < peak.range[1] - half.range) |
                (max(x[peak.interval]) >= max(x)) |
                (min(x[peak.interval]) <= min(x))) {
                boundary <- 1
                break
            } # end if statement
        } # end boundary repeat loop for peak
        if (boundary == 1) {
            peak.time <- NA
            peak <- NA
            peak_loess <- NA
        } else {
            amp.ind <- peak.interval[peak.ind]
            peak.time <- loess.smooth$x[amp.ind] # peak latency
            peak <- mean(y[amp.ind]) # peak amplitude
            peak_loess <- mean(loess.smooth$y[amp.ind]) # loess peak amplitude
        } # end peak boundary if statement
    }# end peak feature
    # Dip feature
    if (find.dip) {
        boundary <- 0
        # Interval for identifying dip
        dip.interval <- which(loess.smooth$x >= dip.range[1] &
                                  loess.smooth$x <= dip.range[2])
        # Check if dip is on the boundary and iterate
        repeat {
            # Identify minimum within interval
            dip.ind <- which.min(loess.smooth$y[dip.interval])
            # Shift interval if maximum on boundary
            if (dip.ind == length(dip.interval)) {
                dip.interval <- dip.interval + 1
            } else if (dip.ind == 1) {
                dip.interval <- dip.interval - 1
            } else {
                break
            } # end boundary if statement for dip
            # Quality check to see if estimate is far from original interval
            half.range <- (dip.range[2] - dip.range[1])/2
            if ((max(x[dip.interval]) > dip.range[2] + half.range) |
                (min(x[dip.interval]) < dip.range[1] - half.range) |
                (max(x[dip.interval]) >= max(x)) |
                (min(x[dip.interval]) <= min(x))) {
                boundary <- 1
                break
            } # end if statement
        } # end boundary repeat loop for dip
        # Identify time
        if (boundary == 1) {
            dip.time <- NA
            dip <- NA
            dip_loess <- NA
        } else {
            amp.ind <- dip.interval[dip.ind]
            dip.time <- loess.smooth$x[amp.ind] # dip latency
            dip <- mean(y[amp.ind]) # dip amplitude
            dip_loess <- mean(loess.smooth$y[dip.ind]) # loess peak amplitude
        } # end dip boundary if statement
    } # end dip feature
    # Output peak,dip magnitude and latency
    if (find.peak & find.dip) {
        features <- c(peak, peak_loess, peak.time, dip,dip_loess, dip.time)
        names(features) <- c("peak", "peak_loess", "peak.time", 
                             "dip", "dip_loess", "dip.time")
    } else if (find.peak & !find.dip) {
        features <- c(peak, peak_loess, peak.time)
        names(features) <- c("peak", "peak_loess", "peak.time")
    } else if (!find.peak & find.dip) {
        features <- c(dip, dip_loess, dip.time)
        names(features) <- c("dip", "dip_loess", "dip.time")
    }
    return(list(loess_smooth = loess_smooth,
                features = features))
}