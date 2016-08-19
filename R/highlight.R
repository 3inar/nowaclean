#' @export
highlight <- function(samplename, pca, box, dens, ma) {
   old.par <- par(mfrow=c(2,2))

   if (!is.character(samplename)) stop("samplename should be character vector")
   if (class(pca) != "prcout") stop("pca class should be prcout")
   if (class(box) != "boxout") stop("box class should be boxout")
   if (class(dens) != "dens") stop("dens class should be dens")
   if (class(ma) != "mapout") stop("ma class should be mapout")

   plot(pca, highlight=samplename)
   plot(box, highlight=samplename)
   plot(dens, highlight=samplename)
   plot(ma, highlight=samplename)

   par(old.par)
}
