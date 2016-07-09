#' @export
mapout <- function(x) {
  obj <- list()
  class(obj) <- "mapout"

  obj$median_array <- plyr::aaply(x, 2, median)

  # m statistics
  obj$M <- plyr::aaply(x, 1, function(array) {
    array - obj$median_array
  })

  # a statistics
  obj$A <- plyr::aaply(x, 1, function(array) {
    (array + obj$median_array)/2
  })

  # Hoeffding
#   obj$hoeffding <- plyr::aaply(1:nrow(x), 1, function(i) {
#       Hmisc::hoeffd(obj$M[i, ], obj$A[i, ])$D[1,2]
#   })

  obj$information <- plyr::aaply(1:nrow(x), 1, function(i) {
    x1 <- obj$A[i,]
    x2 <- obj$M[i,]

    binner <- nclass.FD
    y2d = entropy::discretize2d(x1, x2, numBins1=binner(x1), numBins2=binner(x2))

    # mutual information
    suppressWarnings({
      ret <- entropy::mi.empirical(y2d)
    })

    ret
  })

  obj
}

#' @export
predict.mapout <- function(obj, sdev=2) {
  info <- obj$information
  info > mean(info) + sdev*sd(info)
}

#' @export
plot.mapout <- function(obj, sdev=2) {
  worst <- order(obj$information, decreasing = T)[1]
  plot(obj$A[worst,], obj$M[worst, ])
}
