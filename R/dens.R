#' Estimate array-wise densities
#'
#' This metod basically just runs \code{density()} per row in the input
#' matrix and packages this in an object of class \code{dens}
#'
#' @param x A matrix of array intensities, samples by row.
#'
#' @export
dens <- function(x) {
  if (is.null(rownames(x))) stop("argument x should have rownames")
  ret <- plyr::alply(x, 1, stats::density, .dims=T)
  class(ret) <- "dens"

  ret
}

#' Plot method for \code{dens} objects
#'
#' Plots all the densities in the \code{dens}-object in a single plot.
#'
#' @param obj an object of class "dens"
#' @param batch a vector indicating which batch a sample belongs to. used for coloring by batch.
#' @param highlight sample names of samples to highligth in red on the plot. useful for inspecting specific outliers
#'
#' @export
plot.dens <- function(obj, batch=NULL, highlight=NULL) {
  xrange <- range(plyr::laply(obj, function(o) {o$x}))
  yrange <- range(plyr::laply(obj, function(o) {o$y}))

  plot(NULL, type="n", xlim=xrange, ylim=yrange, ylab="Density", xlab="")

  if (!is.null(highlight)) {
    if (class(highlight) != "character") stop("highlight= should be character")

    out <- obj[names(obj) %in% highlight]
    not <- obj[!names(obj) %in% highlight]
  } else {
    out <- list()
    not <- obj
  }

  if (!is.null(batch)) {
    colv <- RColorBrewer::brewer.pal(12, "Set3")
    batch <- as.numeric(as.factor(batch))
    if (max(batch) > 12) {  # if more than 12 colors, recycle
      batch <- ((batch - 1) %% 12) + 1
    }
    colors <- colv[batch]
    colors <- colors[!names(obj) %in% highlight]
  } else {
    colors <- rep("black", length(not))
  }

  for (i in 1:length(not)) {
    lines(not[[i]], col=colors[i])
  }

  for (d in out) {
    lines(d, col="red")
  }

}

