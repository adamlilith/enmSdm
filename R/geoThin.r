#' Thin geographic points (mostly) deterministically
#'
#' This function thins geographic points such that none have nearest neighbors closer than some user-specified distance. The results are almost deterministic (see Details).
#' @param x Data frame, matrix, SpatialPoints, or SpatialPointsDataFrame object.
#' @param minDist Numeric. Minimum distance needed between points to retain them. Points falling < this distance will be discarded. If \code{distFunct} is \code{distCosine} then this should be in the same units as \code{f} (see \code{link[geosphere]{distCosine}} and related "\code{dist}" functions).
#' @param longLat Two-element character list \emph{or} two-element integer list. If \code{x} is a data frame then this should be a character list specifiying the names of the fields in \code{x} \emph{or} a two-element list of integers that correspond to longitude and latitude (in that order). For example, \code{c('long', 'lat') or \code{c(1, 2)}. If \code{x} is a matrix then this is a two-element list indicating the column numbers in \code{x} that represent longitude and latitude (for example, \code{c(1, 2)}). If \code{x} is a \code{SpatialPoints} or a \code{SpatialPointsDataFrame}  object then this argument is ignored.
#' @param distFunct Either a function or \code{NULL}. If \code{NULL} then \code{\link[geosphere]{distCosine}} is used to calculate distances.  More accurate distances can be obtained by using other functions (see \code{\link[geosphere]{distHaversine}} and references therein). Alternatively, a custom function can be used so long as its first argument is a 2-column numeric matrix with one row for the x- and y-coordinates of a single point and its second argument is a two-column numeric matrix with one or more rows of other points.
#' @param verbose Logical. If \code{TRUE} then display progress.
#' @param na.rm Logical. If TRUE then remove coordinates with \code{NA}s first.
#' @param ... Arguments to pass to \code{distFunct}.
#' @details
#' The procedure for removing points is as follows:
#' \itemize{
#' 	\item Find points with largest number of neighbors (< \code{minDist} away). If just one such point exists, remove it, but if there is more than one then...
#' 	\item Of these find the points with the closest neighbor within \code{minDist}. If just one such point exists, remove it, but if there is more than one then...
#' 	\item Of these find the point that is closest to the centroid of all non-removed points. If just one such point exists, remove it, but if there is more than one...
#' 	\item Of these find the point that has the closest neighbor (even if > \code{minDist}). If just one such point exists, remove it, but if there is more than one then...
#' 	\item Of these find the point that has the smallest median distance to all points (even if > \code{minDist}). If just one such point exists, remove it, but if there is more than one then...
#' 	\item Of these randomly select a point and remove it.
#' 	\item Repeat.
#' }
#' Thus the results are deterministic up to the last tie-breaking step.
#' @return Object of class \code{x}.
#' @seealso \code{\link{geoThinApprox}}
#' @examples
#' x <- data.frame(long=c(-90.1, -90.1, -90.15, -90.17, -90.2, -89),
#'    lat=c(38, 38, 38, 38, 38, 38), point=letters[1:6])
#' geoThin(x, minDist=500, longLat=c(1, 2))
#' geoThin(x, minDist=5000, longLat=c(1, 2))
#' # example of potential randomness
#' set.seed(123)
#' geoThin(x, minDist=1000, longLat=c(1, 2))
#' geoThin(x, minDist=1000, longLat=c(1, 2))
#' geoThin(x, minDist=1000, longLat=c(1, 2))
#' @export

geoThin <- function(
	x,
	minDist,
	longLat = NULL,
	distFunct = NULL,
	verbose = FALSE,
	na.rm = TRUE,
	...
) {

	#################
	### functions ###
	#################

	if (is.null(distFunct)) distFunct <- geosphere::distCosine

	####################
	## pre-processing ##
	####################

	# get coordinates
	xy <- xToCoords(x, longLat)
	index <- seq_along(xy)

	##########
	## MAIN ##
	##########

	# calculate distances between points
	if (verbose) {
		omnibus::say('Calculating pairwise distances...')
		prog <- txtProgressBar(min=0, max=length(xy), width=30, style=3)
	}

	dists <- matrix(rep(NA, length(xy)^2), nrow=length(xy))

	for (i in seq_along(xy)) {

		dists[i, i:ncol(dists)] <- distFunct(xy[i], xy[i:ncol(dists)], ...)
		dists[i, 1:i] <- dists[1:i, i]
		if (verbose) utils::setTxtProgressBar(prog, i)

	}

	if (verbose)omnibus::say('')

	diag(dists) <- Inf

	# remove distances ("dists" will be changed as points are removed)
	masterDists <- dists

	# neighbors
	neighs <- colSums(dists < minDist, na.rm=TRUE)

	# while there is at least one point with neighbors too close
	if (verbose) {
		omnibus::say('Thinning points...')
		totalNeighs <- sum(neighs)
		prog <- utils::txtProgressBar(min=0, max=totalNeighs, width=30, style=3)
	}

	while (any(neighs > 0)) {

		# get points with greatest number of neighbors
		removeThis <- which(neighs == max(neighs))

		# if more than one point has the maximum number of neighbors, break tiena.rm=TRUE
		if (length(removeThis) > 1) {

			# na.rm=TRUE by removing point with closest neighbor
			minDistHighNeigh <- apply(dists[removeThis, removeThis], 2, min, na.rm=TRUE)
			removeThis <- removeThis[which(minDistHighNeigh == min(minDistHighNeigh))]

			# if there is still a tie, break tiena.rm=TRUE
			if (length(removeThis) > 1) {

				# na.rm=TRUE by removing point closest to centroid of all remaining points
				cent <- rgeos::gCentroid(xy[index, ])
				distToCent <- distFunct(cent, xy[removeThis], ...)

				removeThis <- removeThis[which(distToCent == min(distToCent))]

				# if there is *still* a tie, break tiena.rm=TRUE
				if (length(removeThis) > 1) {

					# na.rm=TRUE by removing point with closest nearest neighbor (even if > minDist away)
					closestNeigh <- apply(dists[removeThis, ], 1, min, na.rm=TRUE)
					removeThis <- removeThis[which(closestNeigh == min(closestNeigh))]

					# if there is *still* a tie, break tiena.rm=TRUE
					if (length(removeThis) > 1) {

						# na.rm=TRUE by removing point with least median distance to all other points (even if > minDist away)
						distToAllPoints <- apply(dists[removeThis, ], 1, stats::median, na.rm=TRUE)
						removeThis <- removeThis[which(distToAllPoints == min(distToAllPoints))]

						# if there is *still* a tie, break tiena.rm=TRUE
						if (length(removeThis) > 1) {

							# na.rm=TRUE by randomly selecting a point
							removeThis <- removeThis[sample(removeThis, 1)]

						}

					}

				}

			}

		}

		# remove point
		index <- index[-removeThis]
		dists[removeThis, ] <- rep(Inf, ncol(dists))
		dists[ , removeThis] <- rep(Inf, nrow(dists))

		# neighbors
		neighs <- colSums(dists < minDist, na.rm=TRUE)

		if (verbose) setTxtProgressBar(prog, totalNeighs - sum(neighs))

	}

	if (verbose) omnibus::say('')

	# remove too-neighborly points
	x <- if (class(x) == 'data.frame' | class(x) == 'matrix') {
		x[index, ]
	} else {
		x[index]
	}

	x


}

