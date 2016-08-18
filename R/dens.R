#' @export
dens <- function(x) {
  if (is.null(rownames(x))) stop("argument x should have rownames")
  ret <- plyr::alply(x, 1, density, .dims=T)
  class(ret) <- "dens"

  ret
}

# NB col overrides outliers
#' @export
plot.dens <- function(obj, highlight=NULL) {
  xrange <- range(plyr::laply(obj, function(o) {o$x}))
  yrange <- range(plyr::laply(obj, function(o) {o$y}))

  plot(NULL, type="n", xlim=xrange, ylim=yrange, ylab="Density", xlab="")

  if (!is.null(highlight)) {
    if (class(highlight) != "character") stop("highlight= should be character")

    out <- obj[names(obj) %in% highlight]
    obj <- obj[!names(obj) %in% highlight]

  } else {
    out <- list()
  }
  for (d in obj) {
    lines(d)
  }

  for (d in out) {
    lines(d, col="red")
  }

}

