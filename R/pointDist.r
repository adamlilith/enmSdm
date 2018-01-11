#' Calculate distance matrix between geographic points
#'
#' This function calculates a distance matrix between points. It is similar to the \code{dist()} function except that it does its calculations in geographic space.
#' @param x A data frame, matrix, \code{SpatialPoints}, or \code{SpatialPointsDataFrame} object. If a data frame or matrix then the coordinate reference system is assumed to be unprojected (WGS84).
#' @param distFunct Either a function or \code{NULL}. If \code{NULL} then \code{\link[geosphere]{distCosine}}is used to calculate distances.  More accurate distances can be obtained by using other functions (see \code{\link[geosphere]{distHaversine}} and related "\code{dist}" functions). Alternatively, a custom function can be used so long as its first argument is a 2-column numeric matrix with one row for the x- and y-coordinates of a single point and its second argument is a two-column numeric matrix with one or more rows of other points.
#' @param longLat Two-element character list \emph{or} two-element integer list. If \code{x} is a data frame then this should be a character list specifiying the names of the fields in \code{x} \emph{or} a two-element list of integers that correspond to longitude and latitude (in that order). For example, \code{c('long', 'lat') or \code{c(1, 2)}. If \code{x} is a matrix then this is a two-element list indicating the column numbers in \code{x} that represent longitude and latitude. For example, \code{c(1, 2)}. If \code{x} is a SpatialPoints object then this argument is ignored.
#' @param ... Arguments to pass to \code{distFunct}.
#' @return Matrix of values in units given by the units of \code{f} (if the distance function used is [geosphere::distCosine()] or related \code{dist} functions).
#' @examples
#' \dontrun{
#' set.seed(17)
#' # random points centered on long -90, lat 45
#' x <- cbind(rnorm(5, -90), rnorm(5, 45))
#' pointDist(x)
#' }
#' @export
pointDist <- function(x, distFunct = NULL, longLat = NULL, ...) {

	if (is.null(distFunct)) distFunct <- geosphere::distCosine

	# get coordinates
	xy <- xToCoords(x)

	dists <- matrix(rep(1, length(xy)^2), nrow=length(xy))
	for (i in 1:nrow(dists)) dists[i, i:length(xy)] <- distFunct(xy[i], xy[i:length(xy)], ...)
	dists <- dists + t(upper.tri(dists) * dists)
	dists

}
