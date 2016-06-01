# samples by row
#' @export
boxout <- function(x) {
  obj <- list()
  class(obj) <- "boxout"

  pooledcdf <- ecdf(x)
  statistics <- plyr::aaply(x, 1, function(row) {
    # for plots
    quants <- quantile(row, c(0.25, 0.5, 0.75))
    iqr <- abs(quants[3] - quants[1])
    inner <- abs(row - quants[2]) <= 1.5*iqr
    whisker_l <- min(row[inner])
    whisker_u <- max(row[inner])

    # for outlier detection
    kolmogs <- ks.test(row, pooledcdf)$statistic
    res <- c(whisker_l, quants, whisker_u, kolmogs)
    names(res) <- NULL
    res
  })

  colnames(statistics) <- c("wl", "0.25", "0.5", "0.75", "wu", "ks")

  obj$cdf <- pooledcdf
  obj$statistics <- statistics
  obj
}

#' @export
plot.boxout<- function(x) {
  quant <- x$statistics[, c("wl", "0.25", "0.5", "0.75", "wu")]
  kstats <- x$statistics[, "ks"]
  #quant <- quant[order(quant[, "0.5"]), ]
  quant <- quant[order(kstats), ]
  ymax <- max(quant)
  ymin <- min(quant)

  plot(quant[, "0.5"], type="n", ylim=c(ymin, ymax)) #, ...)

  lines(quant[, "wl"], lty="dashed")
  lines(quant[, "wu"], lty="dashed")
  lines(quant[, "0.25"])
  lines(quant[, "0.75"])
  lines(quant[, "0.5"], lwd=2)
}

#' @export
predict.boxout <- function(obj, sdev=2) {
  kst <- obj$statistics[, "ks"]
  kst - mean(kst) > sdev*sd(kst)
}