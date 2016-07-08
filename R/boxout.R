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
    # i don't care that I can't get exact p-values due to ties, I just
    # want the statistic.
    suppressWarnings(
      kolmogs <- ks.test(row, pooledcdf)$statistic
    )
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
plot.boxout<- function(x, batch=NULL, ...) {
  quant <- x$statistics[, c("wl", "0.25", "0.5", "0.75", "wu")]
  kstats <- x$statistics[, "ks"]

  if (is.null(batch)) {
    quant <- quant[order(kstats), ]
  } else {
    batchorder <- order(batch)
    quant <- quant[batchorder, ]
  }

  ymax <- max(quant)
  ymin <- min(quant)

  plot(quant[, "0.5"], type="n", ylim=c(ymin, ymax)) #, ...)

  if (is.null(batch)) {
    abline(v=(length(kstats) - sum(predict(x, ...))), col="red")
  } else {
    cuts <- which(!duplicated(batch[batchorder]))
    abline(v=cuts[-1], col="gray") # 1st one uninformative
  }

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
