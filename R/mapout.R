#' MA-plot-based outlier detection
#'
#' Computes M and A statistics between each array and a probe-wise median array. If an
#' array is well-behaved, there should be no trend in M as a function of A. This
#' trendiness is estimated by mutual information.
#'
#' @param x a matrix of \strong{log transformed, raw} gene expression values
#'   (not fold change). Observations by row, probes by columns.
#'
#' @return an object of class \code{mapout}
#'
#' @export
mapout <- function(x) {
  if (is.null(rownames(x))) stop("input x should have rownames")

  obj <- list()
  class(obj) <- "mapout"

  obj$median_array <- plyr::aaply(x, 2, stats::median)

  # m statistics
  obj$M <- plyr::aaply(x, 1, function(array) {
    array - obj$median_array
  })

  # a statistics
  obj$A <- plyr::aaply(x, 1, function(array) {
    (array + obj$median_array)/2
  })

  # discretization granularity
  nba <- 200
  nbm <- 200

  ra <- range(obj$A)
  rm <- range(obj$M)

  obj$information <- plyr::aaply(1:nrow(x), 1, function(i) {
    x1 <- obj$A[i,]
    x2 <- obj$M[i,]

    y2d = entropy::discretize2d(x1, x2,
                                numBins1=nba,
                                numBins2=nbm,
                                r1=ra,
                                r2=rm)

    # mutual information
    suppressWarnings({
      ret <- entropy::mi.empirical(y2d)
    })

    ret
  })

  obj
}

#' Predict method for \code{mapout} objects
#'
#' Provides a prediction for whether an observation is a suspected outlier
#'
#' @param object object of class "mapout"
#' @param sdev Number of standard deviations (in mutual information) that an
#'  observation should be from the center to be called an outlier
#' @param ... unused
#'
#' @return a character vector of row names for the outliers as defined by the \code{sdev}
#'  parameter.
#'
#' @export
predict.mapout <- function(object, sdev=3, ...) {
  info <- object$information
  rownames(object$M)[info > mean(info) + sdev*stats::sd(info)]
}

#' Plot method for \code{mapout} objects
#'
#' Plots the MA-plots \code{nout} "worst outliers" in a grid of \code{ncol}
#' columns and plots a loess-smoothed line over the data. You want this line
#' to be as horizontal as possible, but a certain banana-shape is not uncommon.
#' The loess-smoothing might be slow, and you might end up with a lot of points
#' in your plot. To work around this, use \code{subsample=T} (default) to work
#' examine a random subset of 2500 probes.
#'
#' @param x object of class \code{mapout}
#' @param nout now many outliers to show MAplots for
#' @param highlight optional character vector with names of samles to show in the plot.
#'  Overrides the standard plotting of outliers
#' @param lineup show \code{nout} random samples for comparison? defaults to false.
#' @param subsample show random subset of 2500 probes for speed. Defaults to true
#' @param ncol number of columns in output plot
#' @param ... unused
#' @export
plot.mapout <- function(x, nout=1, highlight=NULL, lineup=F, subsample=T,
                        ncol=ifelse(is.null(highlight), nout, length(highlight)),
                        ...) {

  obj <- x
  if (is.null(highlight)) nplots <- ifelse(lineup, 2*nout, nout)
  else nplots <- nplots <- ifelse(lineup, 2*length(highlight), length(highlight))
  nrows <- ceiling(nplots/ncol)

  if (nplots > 1) old.par <- graphics::par(mfrow=c(nrows,ncol))

  if (subsample) {
    nsub <- ceiling(ncol(obj$M))/10
    genes <- sample(1:ncol(obj$M), ifelse(nsub < 2500, nsub, 2500))
  } else {
    genes <- 1:ncol(obj$M)
  }

  yrange <- range(obj$M)
  xrange <- range(obj$A)
  randoms <- sample(rownames(obj$M), nplots/2)

  if (is.null(highlight)) {
    worst <- rownames(obj$M)[order(obj$information, decreasing = T)[1:nout]]

    for (w in worst) {
      single_plot(obj, w, genes, xrange, yrange, spline=T)
    }

  } else {
    for (w in highlight) {
      single_plot(obj, w, genes, xrange, yrange, spline=T)
    }
  }

  if (lineup) {
    for (w in randoms) {
      single_plot(obj, w, genes, xrange, yrange, spline=T, rs=T)
    }
  }

  if (nplots > 1) graphics::par(old.par)
}

fullreport <- function(obj, perpage=20, subsample=T) {
  x <- order(obj$information, decreasing=T)
  pages <-  split(x,ceiling(seq_along(x)/perpage))

  # pdf('eg.pdf', width = 8.3, height = 11.7)  ## Device with dimensions of A4 paper
  # par(mfrow = c(3,2))                        ## 2x3 grid of plotting areas
  # replicate(plot(rnorm(99)), n = 6)          ## 6 example plots
  # dev.off()

  grDevices::pdf(file="report.pdf", width=8.3, height=11.7, onefile=T)
  graphics::par(omi = rep(.5, 4))                      ## 1/2 inch outer margins
  plyr::l_ply(pages, function(candidates) {
    graphics::plot(obj, highlight=candidates, ncol=4, subsample=subsample)
  })
  grDevices::dev.off()
}

outlierplot <- function(obj, sdev=2, ncol=5, subsample=T) {
  outliers <- predict.mapout(obj, sdev=sdev)

  worst <- which(outliers)
  worst <- worst[order(obj$information[outliers], decreasing = T)]

  graphics::plot(obj, highlight=worst, ncol=ncol, subsample=subsample)
}

single_plot <- function(obj, samplename, genes, xrange, yrange, spline=F, rs=F) {
  w <- which(rownames(obj$M)==samplename)

  title <- rownames(obj$M)[w]
  if (rs) title <- paste0(title, " (random)")
  graphics::plot(obj$A[w, genes], obj$M[w, genes], main=title,
       sub=paste0("MI: ", obj$information[w]),
       ylim=yrange, xlim=xrange, xlab="A", ylab="M")
  graphics::abline(h=0, col="grey")

  if (spline) {
    lo <- stats::loess(obj$M[w, genes]~obj$A[w, genes])
    rge <- range(obj$A)
    plotseq <- seq(rge[1], rge[2], length.out=100)
    predictions <- stats::predict(lo, plotseq)
    graphics::lines(plotseq, predictions, col="red")
  }
}
