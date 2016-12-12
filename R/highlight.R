#' Highlight an observation in four different views
#'
#' Provides linking of plots for one observation in the four different
#' models we provide. Puts the plots in a grid an calls \code{plot(x, highlight=observation)}
#' for each of the four models.
#'
#' @param observation the name of the observation to highligh
#' @param pca an object of class "prcout" where the observation occors
#' @param box an object of class "boxout" where the observation occors
#' @param dens an object of class "dens" where the observation occors
#' @param ma an object of class "mapout" where the observation occors
#'
#' @export
highlight <- function(observation, pca, box, dens, ma) {
   old.par <- par(mfrow=c(2,2))

   if (!is.character(observation)) stop("observation should be character vector")
   if (class(pca) != "prcout") stop("pca class should be prcout")
   if (class(box) != "boxout") stop("box class should be boxout")
   if (class(dens) != "dens") stop("dens class should be dens")
   if (class(ma) != "mapout") stop("ma class should be mapout")

   plot(pca, highlight=observation)
   plot(box, highlight=observation)
   plot(dens, highlight=observation)
   plot(ma, highlight=observation)

   par(old.par)
}
