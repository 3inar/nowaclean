#' PCA-based outlier detection
#'
#' Takes a matrix of samples x measurements and looks for outliers in the two
#' first principal components of the data as defined by mahalanobis distance
#' to the center of the data in number of standard deviations
#'
#' @param x a numerical matrix with samples by row, measurements by column
#' @param prob How unlikely should a data point at least be in order to not be
#'  considered part of the "center mass" of the data. Translated to k in
#'  Chebyshev's inequality P(|Z| >= k) =< 1/k^2 and applied to the two
#'  first PCs.
#'
#' @return an object of class \code{prcout}
#'
#' @export
prcout <- function(x, prob=0.01) {
  if (is.null(rownames(x))) stop("input x should have rownames")
  obj <- list()
  obj$prc <- stats::prcomp(x)

  # identify inner 100(1-prob)% mass of data in PCs 1 and 2
  prc <- obj$prc$x[, 1:2]
  alpha <- prob
  k <- sqrt(1/alpha)
  x95 <- abs((prc[, 1] - mean(prc[, 1]))/stats::sd(prc[, 1])) < k
  y95 <- abs((prc[, 2] - mean(prc[, 2]))/stats::sd(prc[, 2])) < k
  prc_inner <- prc[x95 & y95, ]

  # define mahalanobis distances af f'n of inner data, apply to all data
  obj$mu <- colMeans(prc_inner)
  obj$Sigma <- stats::cov(prc_inner)
  obj$sd <- stats::sd(stats::mahalanobis(prc_inner, obj$mu, obj$Sigma))
  obj$mahalanobis <- stats::mahalanobis(prc, obj$mu, obj$Sigma)

  class(obj) <- "prcout"
  obj
}

#' Plot method for \code{prcout} objects
#'
#' Plots the two first PCs of the data with contour lines corresponding
#' to the mahalanobis distance to the data center (as defined by the inner mass
#' of the data). Highlights outliers or the samples defined in \code{highlight=}.
#' Alternately colors points by batch if the \code{batch} parameter is defined
#'
#' @param x object of class \code{prcout}
#' @param batch optional vector that indicates which bach a sample belogs to.
#'  Points will be colored by batch if this vector is provided. Overrides
#'  \code{highlight=}
#' @param highlight optional character vector with names of samles to highlight in the plot.
#'  Overrides the highlighting of outliers
#' @param ... passed to predict(obj, ...) to label outliers
#' @seealso \code{\link{prcout}}, \code{\link{predict.prcout}}
#' @export
plot.prcout <- function(x, batch=NULL, highlight=NULL, ...) {
  obj <- x
  prc <- obj$prc$x[, 1:2]

  # grid for contour lines
  r1 <- range(prc[,1])
  r2 <- range(prc[,2])
  s1 <- seq(r1[1], r1[2], length.out=50)
  s2 <- seq(r2[1], r2[2], length.out=50)
  grid <- expand.grid(s1,s2)
  mha <- matrix(stats::mahalanobis(grid, obj$mu, obj$Sigma), nrow=length(s1))

  # plot
  graphics::contour(s1, s2, mha/obj$sd, levels=1:6)

  if (is.null(highlight)) {
    hset <- predict.prcout(obj, ...)
  }  else {
    hset <- highlight
  }

  if (!is.null(batch)) {
    colv <- RColorBrewer::brewer.pal(12, "Set3")
    batch <- as.numeric(as.factor(batch))
    if (max(batch) > 12) {  # if more than 12 colors, recycle
      batch <- ((batch - 1) %% 12) + 1
    }
    colors <- colv[batch]
    graphics::points(prc, col=colors, pch=20)
  }  else {
    colors <- rep("black", nrow(prc))
    colors[rownames(prc) %in% hset] <- "red"
    graphics::points(prc[order(colors), ], col=colors[order(colors)], pch=20)
  }

}

#' Predict method for \code{prcout} objects
#'
#' Provides a prediction for whether an observation is a suspected outlier
#'
#' @param object object of class "prcout"
#' @param sdev Number of standard deviations (in mahalanobis distance) that an
#'  observation should be from the center to be called an outlier
#' @param ... unused
#'
#' @return a character vector of row names for the outliers as defined by the \code{sdev}
#'  parameter.
#'
#' @export
predict.prcout <- function(object, sdev=3, ...) {
  rownames(object$prc$x)[object$mahalanobis > sdev*object$sd]
}

