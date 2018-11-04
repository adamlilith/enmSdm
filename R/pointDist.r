#' Calculate distance matrix between geographic points
#'
#' This function calculates a distance matrix between points. It is similar to the \code{dist()} function except that it does its calculations in geographic space.
#' @param x A data frame, matrix, \code{SpatialPoints}, or \code{SpatialPointsDataFrame} object. If a data frame or matrix then the coordinate reference system is assumed to be unprojected (WGS84).
#' @param y If \code{NULL} then pairwise distances are calculated between all points in \code{x}. However, if \code{y} is a data frame, matrix, \code{SpatialPoints}, or \code{SpatialPointsDataFrame} object then distances are between all points in \code{x} and \code{y}.
#' @param distFunct Either a function or \code{NULL}. If \code{NULL} then \code{\link[geosphere]{distCosine}}is used to calculate distances.  More accurate distances can be obtained by using other functions (see \code{\link[geosphere]{distHaversine}} and related "\code{dist}" functions). Alternatively, a custom function can be used so long as its first argument is a 2-column numeric matrix with one row for the x- and y-coordinates of a single point and its second argument is a two-column numeric matrix with one or more rows of other points.
#' @param longLat Two-element character list \emph{or} two-element integer list. If \code{x} (and \code{y}) is a data frame then this should be a character list specifying the names of the fields in \code{x} \emph{or} a two-element list of integers that correspond to longitude and latitude (in that order). For example, \code{c('long', 'lat')} or \code{c(1, 2)}. If \code{x} (and \code{y}) is a matrix then this is a two-element list indicating the column numbers in \code{x} (and \code{y}) that represent longitude and latitude. For example, \code{c(1, 2)}. If \code{x} is a \code{SpatialPoints} object then this argument is ignored.
#' @param ... Arguments to pass to \code{distFunct}.
#' @return Matrix of values in units given by the coordinate reference system of \code{x} or in meters (if \code{x} is a matrix, data frame, or unprojected).
#' @seealso \code{\link[geosphere]{distCosine}}, \code{\link[omnibus]{pairDist}}
#' @examples
#' set.seed(17)
#' # random points centered on long -90, lat 45
#' x <- cbind(rnorm(5, -90), rnorm(5, 45))
#' pointDist(x)
#' @export
pointDist <- function(x, y = NULL, distFunct = NULL, longLat = NULL, ...) {

	if (is.null(distFunct)) distFunct <- geosphere::distCosine

	# get coordinates
	x <- xToCoords(x, longLat=longLat)
	if (is.null(y)) {
		y <- x
	} else {
		y <- xToCoords(y, longLat=longLat)
	}

	dists <- matrix(NA, nrow=length(x), ncol=length(y))
	for (i in 1:nrow(dists)) dists[i, ] <- distFunct(x[i], y, ...)

	dists

}
