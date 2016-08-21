#' @export
dens <- function(x) {
  if (is.null(rownames(x))) stop("argument x should have rownames")
  ret <- plyr::alply(x, 1, density, .dims=T)
  class(ret) <- "dens"

  ret
}

# NB col overrides outliers
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

