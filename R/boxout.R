#' Outlier detection and batch-effect investigation by boxplots
#'
#' Compared the empirical distribution function of each array's intensities
#' to the (smoothed) ecdf of all arrays' intensities pooled. The comparison
#' is done by the Kolmogorov--Smirnov statistic where we ignore the presence of
#' ties.
#'
#' @param x A matrix of intensities, observations by row.
#'
#' @return An object of class "boxout"
#' @export
boxout <- function(x) {
  if (is.null(rownames(x))) stop("input x should have rownames")

  obj <- list()
  class(obj) <- "boxout"

  pooledcdf <- stats::ecdf(x)
  statistics <- plyr::aaply(x, 1, function(row) {
    # for plots
    quants <- stats::quantile(row, c(0.25, 0.5, 0.75))
    iqr <- abs(quants[3] - quants[1])
    inner <- abs(row - quants[2]) <= 1.5*iqr
    whisker_l <- min(row[inner])
    whisker_u <- max(row[inner])

    # for outlier detection
    # i don't care that I can't get exact p-values due to ties, I just
    # want the statistic.
    suppressWarnings(
      kolmogs <- stats::ks.test(row, pooledcdf)$statistic
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
#' (see \code{boxout} documentation) and shows a red line that marks the
#' first outlier (see \code{predict.boxout} documenation).
#'
#' @param x an object of class \code{boxout}
#' @param batch optional vector of batch labels for each observation. If provided,
#'  the observations will be sorted by batch and the batches separated by grey lines.
#' @param highlight optional character vector of observations to highlight with a
#'    red line. Overrides the highlighting of first outlier.
#' @param ... passed to predict(obj, ...) to label outliers
#' @seealso \code{\link{boxout}}, \code{\link{predict.boxout}}
#' @export
plot.boxout<- function(x, batch=NULL, highlight=NULL, ...) {
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

  graphics::plot(quant[, "0.5"], type="n", ylim=c(ymin, ymax), xlab="Arrays",
                 ylab="Intensity", xaxt="n") #, ...)

  if (!is.null(highlight)) {
    hli <- which(rownames(quant) %in% highlight)
  } else if(is.null(batch)) {
    hli <- (length(kstats) - length(predict.boxout(x, ...)))
  } else {
    hli <- c()
  }

  if (!is.null(batch)) {
    cuts <- which(!duplicated(batch[batchorder]))
    graphics::abline(v=cuts[-1], col="gray") # 1st one uninformative
  }

  graphics::abline(v=hli, col="red")

  graphics::lines(quant[, "wl"], lty="dashed")
  graphics::lines(quant[, "wu"], lty="dashed")
  graphics::lines(quant[, "0.25"])
  graphics::lines(quant[, "0.75"])
  graphics::lines(quant[, "0.5"], lwd=2)

}

#' Predict method for \code{boxout} objects
#'
#' Provides a prediction for whether an observation is a suspected outlier
#'
#' @param object an object of class "boxout"
#' @param sdev Number of standard deviations (in KS statistic) that an
#'  observation should be larger than the mean KS statistic in order to be
#'  considered an outlier
#' @param ... unused
#'
#' @return a vector of names for the outlying observations as defined by the \code{sdev}
#'  parameter.
#'
#' @export
predict.boxout <- function(object, sdev=3, ...) {
  kst <- object$statistics[, "ks"]
  rownames(object$statistics)[kst - mean(kst) > sdev*stats::sd(kst)]
}
