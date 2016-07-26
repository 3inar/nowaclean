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
plot.mapout <- function(obj, sdev=2, nout=1, lineup=F, subsample=T) {
  nrows <- ifelse(lineup, 2, 1)
  old.par <- par(mfrow=c(nrows,nout))

  if (subsample) {
    genes <- sample(1:ncol(obj$M), 2500)
  } else {
    genes <- 1:ncol(obj$M)
  }

  worst <- order(obj$information, decreasing = T)[1:nout]
  randoms <- sample(1:nrow(obj$M), nout)

  if (lineup) {
    yrange <- range(obj$M[unique(c(worst, randoms)), ])
  } else {
    yrange <- range(obj$M[worst, ])
  }

  for (w in worst) {
    plot(obj$A[w, genes], obj$M[w, genes], main=rownames(obj$M)[w], ylim=yrange)
    abline(h=0, col="red")
  }

  if (lineup) {
    for (w in randoms) {
      plot(obj$A[w, genes], obj$M[w, genes], main="random sample", ylim=yrange)
      abline(h=0, col="red")
    }
  }

  par(old.par)
}
