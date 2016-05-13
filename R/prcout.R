#' PCA-based outlier detection
#'
#' Takes a matrix of samples x measurements and looks for outliers in the two
#' first principal components of the data as defined by mahalanobis distance
#' to the center of the data in number of standard deviations
#'
#' @param x a numerical matrix with samples by row, measurements by column
#' @param prob How unlikely should a data point at least be in order to not be
#'  considered part of the "center mass" of the data. Translated to $k$ in
#'  Chebyshev's inequality $P(|Z| \geq k) \leq 1/k^2$ and applied to the two
#'  first PCs.
#'
#' @return an object of class \code{prcout}. Can be plottet with plot(obj)
#'
#' @export
prcout <- function(x, prob=0.01) {
  obj <- list()
  obj$prc <- prcomp(x)

  # identify inner 100(1-prob)% mass of data in PCs 1 and 2
  prc <- obj$prc$x[, 1:2]
  alpha <- prob
  k <- sqrt(1/alpha)
  x95 <- abs((prc[, 1] - mean(prc[, 1]))/sd(prc[, 1])) < k
  y95 <- abs((prc[, 2] - mean(prc[, 2]))/sd(prc[, 2])) < k
  prc_inner <- prc[x95 & y95, ]

  # define mahalanobis distances af f'n of inner data, apply to all data
  obj$mu <- colMeans(prc_inner)
  obj$Sigma <- cov(prc_inner)
  obj$sd <- sd(mahalanobis(prc_inner, obj$mu, obj$Sigma))
  obj$mahalanobis <- mahalanobis(prc, obj$mu, obj$Sigma)

  class(obj) <- "prcout"
  obj
}

#' Plot method for \code{prcout} objects
#'
#' Plots the two first PCs of the data with contour lines corresponding
#' to the mahalanobis distance to the data center (as defined by the inner mass
#' of the data)
#'
#' @param ... passed to predict(obj, ...) to label outliers
#' @seealso \code{\link{prcout}}, \code{\link{predict.prcout}}
#' @export
plot.prcout <- function(obj, ...) {
  prc <- obj$prc$x[, 1:2]

  # grid for contour lines
  r1 <- range(prc[,1])
  r2 <- range(prc[,2])
  s1 <- seq(r1[1], r1[2], length.out=50)
  s2 <- seq(r2[1], r2[2], length.out=50)
  grid <- expand.grid(s1,s2)
  mha <- matrix(mahalanobis(grid, obj$mu, obj$Sigma), nrow=length(s1))

  # plot
  contour(s1, s2, mha/obj$sd, levels=1:6)
  points(prc, col=ifelse(predict(obj, ...), "red", "black"))
}

#' Predict outliers
#'
#' Provides a prediction for whether an observation is a suspected outlier
#'
#' @param sdev Number of standard deviations (in mahalanobis distance) that an
#'  observation should be from the center to be called an outlier
#'
#' @return a logical vector indicating outlier Y/N as defined by the \code{sdev}
#'  parameter.
#'
#' @export
predict.prcout <- function(obj, sdev=2) {
  obj$outliers <- obj$mahalanobis > sdev*obj$sd
}

