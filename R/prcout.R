#' PCA-based outlier detection
#'
#' Takes a matrix of samples x measurements and looks for outliers in the two
#' first principal components of the data as defined by mahalanobis distance
#' to the center of the data in number of standard deviations
#'
#' @param x a numerical matrix with samples by row, measurements by column
#' @param sdev Number of standard deviations (in mahalanobis distance) that an
#'  observation should be from the center to be called an outlier
#'
#' @return an object of class \code{prcout}. Can be plottet with plot(obj)
#'
#' @export
prcout <- function(x, sdev=2) {
  obj <- list()
  obj$prc <- prcomp(x)

  # mahalanobis distances
  prc <- obj$prc$x[, 1:2]
  mu <- colMeans(prc)
  Sigma <- cov(prc)
  obj$mahalanobis <- mhd <- mahalanobis(prc, mu, Sigma)
  sdev_distance <- sd(mhd)
  obj$outliers <- mhd > sdev*sdev_distance

  class(obj) <- "prcout"
  obj
}

#' @export
plot.prcout <- function(obj) {
  prc <- obj$prc$x[, 1:2]
  mu <- colMeans(prc)
  Sigma <- cov(prc)
  r1 <- range(prc[,1])
  r2 <- range(prc[,2])
  s1 <- seq(r1[1], r1[2], length.out=50)
  s2 <- seq(r2[1], r2[2], length.out=50)
  grid <- expand.grid(s1,s2)
  mha <- matrix(mahalanobis(grid, mu, Sigma), nrow=length(s1))
  sdev_distance <- sd(obj$mahalanobis)
  contour(s1, s2, mha/sdev_distance, levels=1:6)
  points(prc, col=ifelse(obj$outliers, "red", "black"))
}

