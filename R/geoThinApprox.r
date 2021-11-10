#' Thin geographic points stochastically
#'
#' This function thins a set of geographic points so of the remainder, none are closer than a given distance. The function is nearly the same as the function \code{thin.algorithm} in the \pkg{spThin} package, except that it accepts a data frame, matrix, SpatialPoints, or SpatialPointsDataFrame as a main argument and the user can specify the distance function to be used and that it returns an object of the same kind.
#' @param x Data frame, matrix, or SpatialPoints* object. See \code{Details} for further information on the coordinate reference system.
#' @param minDist Numeric. Minimum distance (usually in m) thinned points must be from their nearest neighbor.
#' @param longLat Two-element character list \emph{or} two-element integer list. If \code{x} is a data frame then this should be a character list specifiying the names of the fields in \code{x} \emph{or} a two-element list of integers that correspond to longitude and latitude (in that order). For example, \code{c('long', 'lat')} or \code{c(1, 2)}. If \code{x} is a matrix then this is a two-element list indicating the column numbers in \code{x} that represent longitude and latitude. For example, \code{c(1, 2)}. If \code{x} is a \code{SpatialPoints} object then this argument is ignored.
#' @param distFunct Either a function or \code{NULL}. If \code{NULL} then \code{\link[geosphere]{distGeo}} is used to calculate distances.  More accurate distances can be obtained by using other functions (see \code{\link[geosphere]{distHaversine}} and references therein). Alternatively, a custum function can be used so long as its first argument is a 2-column numeric matrix with one row for the x- and y-coordinates of a single point and its second argument is a two-column numeric matrix with one or more rows of other points.
#' @param verbose Logical. If \code{TRUE} then display progress.
#' @param ... Extra arguments to pass to \code{distFunct}.
#' @return Object of class \code{x}.
#' @details If \code{x} is a data frame or a matrix then it will be assumed to be unprojected (WGS84 coordinate reference system) and \code{minDist} should be in units of meters unless the argument \code{r} (passed to the distance function using \code{...}, see \code{\link[geosphere]{distGeo}} is not in meters.
#' @seealso \code{\link[geosphere]{distGeo}}, \code{\link{geoThin}}
#' @examples
#' x <- data.frame(long=c(-90.1, -90.1, -90.15, -90.17, -90.2, -89),
#'    lat=c(38, 38, 38, 38, 38, 38), point=letters[1:6])
#' set.seed(123)
#' geoThinApprox(x, 10000, longLat=c(1, 2)) # run #1
#' geoThinApprox(x, 10000, longLat=c(1, 2)) # run #2
#' geoThinApprox(x, 10000, longLat=c(1, 2)) # run #3
#' geoThinApprox(x, 10, longLat=c(1, 2))
#'
#' # example using SpatialPointsDataFrame
#' data(lemur)
#' data(mad0)
#' par(mfrow=c(1, 3))
#' for (count in 1:3) {
#'     plot(mad0, main='Madagascar')
#'     points(lemur, col='red')
#'     thinned <- geoThinApprox(lemur, 50000)
#'     points(thinned, pch=16)
#'     legend('topright', legend=c('retained', 'discarded'),
#'     col=c('black', 'red'), pch=c(16, 1))
#' }
#' @export

geoThinApprox <- function(
	x,
	minDist,
	longLat = NULL,
	distFunct = NULL,
	verbose = FALSE,
	...
) {

	if (is.null(distFunct)) distFunct <- geosphere::distGeo

	# get coordinates
	xy <- xToCoords(x, longLat)
	index <- seq_along(xy)

	if (verbose) omnibus::say('Calculating pairwise distances between points...')

	# identify points with neighbors within minimum distance
	dists <- geosphere::distm(xy, fun=distFunct)

	diag(dists) <- NA
	isNeigh <- dists <= minDist

	# vector to store number of points too close to each particular point
	neighs <- colSums(isNeigh, na.rm=TRUE)
	neighsStart <- sum(neighs > 0)
	
	# nothing to thin!
	if (neighsStart == 0) {
	
		if (verbose) omnibus::say('No points are <= ', minDist, ' from one another. 
			 Thinning not performed by function geoThinApprox.')
	
	# thin!
	} else {
		
		if (verbose) {
			omnibus::say('Thinning points...')
			pointsWithNeighs <- sum(neighsStart > 0)
			prog <- utils::txtProgressBar(min=0, max=neighsStart, width=32, style=3)
		}

		remainingPointsWithNeighs <- sum(neighs > 0)

		while (sum(!is.na(index)) > 1 & remainingPointsWithNeighs > 0) {

			badNeighs <- which(neighs > 0)
			badIndex <- seq_along(badNeighs)
			remove <- badNeighs[sample(badIndex, 1)]

			index[remove] <- NA

			isNeigh[remove, ] <- FALSE
			isNeigh[ , remove] <- FALSE

			neighs <- colSums(isNeigh, na.rm=TRUE)

			remainingPointsWithNeighs <- sum(neighs > 0)
			
			if (verbose) utils::setTxtProgressBar(prog, neighsStart - remainingPointsWithNeighs)

		}

		index <- index[!is.na(index)]
		
	} # if thinning
		
	# remove too-neighborly points
	x <- if (class(x) == 'data.frame' | class(x) == 'matrix') {
		x[index, ]
	} else if (class(x) == 'SpatialPointsDataFrame') {
		x[index, ]
	} else {
		x[index]
	}

	x

}
