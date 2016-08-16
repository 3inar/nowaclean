#' @export
dens <- function(expressions) {
  ret <- plyr::alply(expressions, 1, density, .dims=T)
  class(ret) <- "dens"

  ret
}

# NB col overrides outliers
#' @export
plot.dens <- function(obj, outliers=NULL) {
  xrange <- range(plyr::laply(obj, function(o) {o$x}))
  yrange <- range(plyr::laply(obj, function(o) {o$y}))

  plot(NULL, type="n", xlim=xrange, ylim=yrange, ylab="Density", xlab="")

  if (!is.null(outliers)) {
    out <- obj[names(obj) %in% out]
    obj <- obj[!names(obj) %in% out]

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

