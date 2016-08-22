# NB assumes log transformed
#' @export
mapout <- function(x) {
  if (is.null(rownames(x))) stop("input x should have rownames")

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

#' @export
predict.mapout <- function(obj, sdev=3) {
  info <- obj$information
  rownames(obj$M)[info > mean(info) + sdev*sd(info)]
}

#' @export
plot.mapout <- function(obj, nout=1, highlight=NULL, lineup=F, subsample=T,
                        ncol=ifelse(is.null(highlight), nout, length(highlight))) {

  if (is.null(highlight)) nplots <- ifelse(lineup, 2*nout, nout)
  else nplots <- nplots <- ifelse(lineup, 2*length(highlight), length(highlight))
  nrows <- ceiling(nplots/ncol)

  if (nplots > 1) old.par <- par(mfrow=c(nrows,ncol))

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

  if (nplots > 1) par(old.par)
}

fullreport <- function(obj, perpage=20, subsample=T) {
  x <- order(obj$information, decreasing=T)
  pages <-  split(x,ceiling(seq_along(x)/perpage))

  # pdf('eg.pdf', width = 8.3, height = 11.7)  ## Device with dimensions of A4 paper
  # par(mfrow = c(3,2))                        ## 2x3 grid of plotting areas
  # replicate(plot(rnorm(99)), n = 6)          ## 6 example plots
  # dev.off()

  pdf(file="report.pdf", width=8.3, height=11.7, onefile=T)
  par(omi = rep(.5, 4))                      ## 1/2 inch outer margins
  plyr::l_ply(pages, function(candidates) {
    plot(obj, highlight=candidates, ncol=4, subsample=subsample)
  })
  dev.off()
}

outlierplot <- function(obj, sdev=2, ncol=5, subsample=T) {
  outliers <- predict(obj, sdev=sdev)

  worst <- which(outliers)
  worst <- worst[order(obj$information[outliers], decreasing = T)]

  plot(obj, highlight=worst, ncol=ncol, subsample=subsample)
}

single_plot <- function(obj, samplename, genes, xrange, yrange, spline=F, rs=F) {
  w <- which(rownames(obj$M)==samplename)

  title <- rownames(obj$M)[w]
  if (rs) title <- paste0(title, " (random)")
  plot(obj$A[w, genes], obj$M[w, genes], main=title,
       sub=paste0("MI: ", obj$information[w]),
       ylim=yrange, xlim=xrange)
  abline(h=0, col="red")

  if (spline) {
    lo <- loess(obj$M[w, genes]~obj$A[w, genes])
    curve(predict(lo, x), col="red", lty=2, add=T)
  }
}
