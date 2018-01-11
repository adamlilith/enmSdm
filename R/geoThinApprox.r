#' Thins a geographic points such that none are within a minimum distance of their nearest neighbor
#'
#' This function thins a set of geographic points so of the remainder, none are closer than a given distance. The function is nearly the same as \code{thin.algorithm()} in the \code{spThin} package, except that it accepts a data frame, matrix, SpatialPoints, or SpatialPointsDataFrame as a main argument and the user can specify the distance function to be used.  Its key advantage over \code{thin.algorithm()} is that 1) it returns the points plus any associated data, whereas that function only returns points; and 2) it is faster, especially for large data sets.
#' @param x Data frame, matrix, or SpatialPoints* object. See \code{Details} for further information on the coordinate reference system.
#' @param minDist Numeric. Minimum distance (usually in m) thinned points must be from their nearest neighbor.
#' @param longLat Two-element character list \emph{or} two-element integer list. If \code{x} is a data frame then this should be a character list specifiying the names of the fields in \code{x} \emph{or} a two-element list of integers that correspond to longitide and latitude (in that order). For example, \code{c('long', 'lat') or \code{c(1, 2)}. If \code{x} is a matrix then this is a two-element list indicating the column numbers in \code{x} that repredsent longitude and latitude. For example, \code{c(1, 2)}. If \code{x} is a \code{SpatialPoints} object then this argument is ignored.
#' @param distFunct Either a function or \code{NULL}. If \code{NULL} then \code{\link[geosphere]{distCosine}} is used to calculate distances.  More accurate distances can be obtained by using other functions (see \code{\link[geosphere]{distHaversine}} and references therein). Alternatively, a custum function can be used so long as its first argument is a 2-column numeric matrix with one row for the x- and y-coordinates of a single point and its second argument is a two-column numeric matrix with one or more rows of other points.
#' @param verbose Logical. If \code{TRUE} then display progress.
#' @param ... Extra arguments to pass to \code{distFunct}.
#' @return Object of class \code{x}.
#' @details If \code{x} is a data frame or a matrix then it will be assumed to be unprojected (WGS84 coordinate reference system) and \code{minDist} should be in units of meters unless the argument \code{r} (passed to the distance function using \code{...}, see \code{\link[geosphere]{distCosine}} is not in meters.
#' @seealso \code{\link[ENMEval]{thin.algorithm}}, \code{\link[geosphere]{distCosine}}, \code{\link{geoThin}}
#' @examples
#' x <- data.frame(long=c(-90.1, -90.1, -90.2, 20), lat=c(38, 38, 38, 38), point=letters[1:4])
#' set.seed(123)
#' geoThinApprox(x, 10000, longLat=c(1, 2)) # run #1
#' geoThinApprox(x, 10000, longLat=c(1, 2)) # run #2
#' geoThinApprox(x, 10000, longLat=c(1, 2)) # run #3
#' geoThinApprox(x, 10000, longLat=c(1, 2)) # run #3
#' geoThinApprox(x, 10, longLat=c(1, 2))
#' @export

geoThinApprox <- function(
	x,
	minDist,
	longLat = NULL,
	distFunct = NULL,
	verbose = FALSE,
	...
) {

#################
### functions ###
#################

if (is.null(distFunct)) distFunct <- geosphere::distCosine

######################
### pre-processing ###
######################

	# get coordinates
	xy <- xToCoords(x, longLat)
	index <- seq_along(xy)

	# calculate distances between points
	dists <- matrix(rep(NA, length(xy)^2), nrow=length(xy))

	if (verbose) {
		omnibus::say('Calculating pairwise distances...')
		prog <- utils::txtProgressBar(min=0, max=length(xy), width=30, style=3)
	}

	for (i in seq_along(xy)) {

		dists[i, i:ncol(dists)] <- distFunct(xy[i], xy[i:ncol(dists)], ...)
		dists[i, 1:i] <- dists[1:i, i]
		if (verbose) utils::setTxtProgressBar(prog, i)

	}

	diag(dists) <- NA
	dists <- dists < minDist

	###################
	### thin points ###
	###################

	# vector to store number of points too close to each particular point
	tooClose <- tooCloseStart <- colSums(dists, na.rm=T)

	if (verbose) {
		omnibus::say('Thinning points...')
		prog <- txtProgressBar(min=0, max=tooCloseStart, width=30, style=3)
	}

	while (length(index) > 1 & sum(tooClose) > 0) {

		# indices of points too close to oner another
		indexMaxNeigh <- index[tooClose == max(tooClose)]

		# if just one point too close to others, sample among them
		if (length(indexMaxNeigh) > 1) indexMaxNeigh <- indexMaxNeigh[sample(length(indexMaxNeigh), 1)]

		# remove point
		index <- index[-which(index==indexMaxNeigh)]
		tooClose <- colSums(dists[index, index], na.rm=T)

		if (verbose) setTxtProgressBar(prog, tooCloseStart - tooClose)

	}


	# remove too-neighborly points
	x <- if (class(x) == 'data.frame' | class(x) == 'matrix') {
		x[index, ]
	} else {
		x[index]
	}

	x

}
