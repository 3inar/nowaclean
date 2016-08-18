#' Outlier detection and batch-effect investigation by boxplots
#'
#' Compared the empirical distribution function of each array's intensities
#' to the (smoothed) ecdf of all arrays' intensities pooled. The comparison
#' is done by the Kolmogorov--Smirnov statistic where we ignore the presence of
#' ties.
#'
#' @param x A matrix of array intensities, samples by row.
#'
#' @return An object of class "boxout"
#' @export
boxout <- function(x) {
  if (is.null(rownames(x))) stop("input x should have rownames")

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

#' Plot method for \code{boxout} objects
#'
#' Plot Karl Broman-style boxplots. By default orders them by KS-statistic
#' (see \code{boxout} documentation) and shows a red line that indicates
#' outlying arrays (see \code{predict.boxout} documenation).
#'
#' @param x an object of class \code{boxout}
#' @param batch optional vector of batch labels for each array. If provided,
#'  the arrays will be sorted by batch and the batches separated by grey lines.
#' @param ... passed to predict(obj, ...) to label outliers
#' @seealso \code{\link{boxout}}, \code{\link{predict.boxout}}
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

#' Predict method for \code{boxout} objects
#'
#' Provides a prediction for whether an observation is a suspected outlier
#'
#' @param sdev Number of standard deviations (in KS statistic) that an
#'  observation should be larger than the mean KS statistic in order to be
#'  considered an outlier
#'
#' @return a logical vector indicating outlier Y/N as defined by the \code{sdev}
#'  parameter.
#'
#' @export
predict.boxout <- function(obj, sdev=3) {
  kst <- obj$statistics[, "ks"]
  kst - mean(kst) > sdev*sd(kst)
}
