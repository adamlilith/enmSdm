#' Thin geographic points (mostly) deterministically
#'
#' This function thins geographic points such that none have nearest neighbors closer than some user-specified distance. The results are almost deterministic (see Details).
#' @param x Data frame, matrix, SpatialPoints, or SpatialPointsDataFrame object.
#' @param minDist Numeric. Minimum distance needed between points to retain them. Points falling < this distance will be discarded. If \code{distFunct} is \code{distGeo} then this should be in the same units as \code{f} (see \code{link[geosphere]{distGeo}} and related "\code{dist}" functions).
#' @param longLat Two-element character list \emph{or} two-element integer list. If \code{x} is a data frame then this should be a character list specifiying the names of the fields in \code{x} \emph{or} a two-element list of integers that correspond to longitude and latitude (in that order). For example, \code{c('long', 'lat')} or \code{c(1, 2)}. If \code{x} is a matrix then this is a two-element list indicating the column numbers in \code{x} that represent longitude and latitude (for example, \code{c(1, 2)}). If \code{x} is a \code{SpatialPoints} or a \code{SpatialPointsDataFrame}  object then this argument is ignored.
#' @param distFunct Either a function or \code{NULL}. If \code{NULL} then \code{\link[geosphere]{distGeo}} is used to calculate distances.  More accurate distances can be obtained by using other functions (see \code{\link[geosphere]{distHaversine}} and references therein). Alternatively, a custom function can be used so long as its first argument is a 2-column numeric matrix with one row for the x- and y-coordinates of a single point and its second argument is a two-column numeric matrix with one or more rows of other points.
#' @param verbose Logical. If \code{TRUE} then display progress.
#' @param ... Arguments to pass to \code{distFunct}.
#' @details
#' The procedure for removing points is as follows:
#' \itemize{
#' 	\item Find points with largest number of neighbors (< \code{minDist} away). If just one such point exists, remove it, but if there is more than one then...
#' 	\item Of these find the points with the closest neighbor within \code{minDist}. If just one such point exists, remove it, but if there is more than one then...
#' 	\item Of these find the point that is closest to the centroid of all non-removed points. If just one such point exists, remove it, but if there is more than one...
#' 	\item Of these find the point that has the smallest median distance to all points (even if > \code{minDist}). If just one such point exists, remove it, but if there is more than one then...
#' 	\item Of these randomly select a point and remove it.
#' 	\item Repeat.
#' }
#' Thus the results are deterministic up to the last tie-breaking step.
#' @return Object of class \code{x}.
#' @seealso \code{\link{geoThinApprox}}
#' @examples
#' # example using data frame
#' x <- data.frame(long=c(-90.1, -90.1, -90.15, -90.17, -90.2, -89),
#'    lat=c(38, 38, 38, 38, 38, 38), point=letters[1:6])
#' x
#' geoThin(x, minDist=500, longLat=1:2, verbose=TRUE)
#' geoThin(x, minDist=5000, longLat=c(1, 2), verbose=TRUE)
#'
#' # example of potential randomness
#' set.seed(111)
#' geoThin(x, minDist=1000, longLat=c(1, 2))
#' geoThin(x, minDist=1000, longLat=c(1, 2))
#' geoThin(x, minDist=1000, longLat=c(1, 2))
#'
#' # example using SpatialPointsDataFrame
#' data(lemurs)
#' fulvus <- lemurs[lemurs$species == 'Eulemur fulvus', c('longitude', 'latitude')]
#' fulvus <- sp::SpatialPointsDataFrame(
#' 		fulvus,
#' 		data=fulvus,
#' 		proj4string=getCRS('wgs84', TRUE)
#' )
#'
#' data(mad0)
#' sp::plot(mad0, main='Madagascar')
#' points(fulvus, col='red')
#' thinned <- geoThin(fulvus, 50000)
#' points(thinned, pch=16)
#' legend('topright', legend=c('retained', 'discarded'),
#' col=c('black', 'red'), pch=c(16, 1))
#'
#' # test to see function works when no points need removed
#' thinned <- geoThin(fulvus, 200, verbose=TRUE)
#' sp::plot(mad0, main='Madagascar')
#' points(fulvus, col='red')
#' points(thinned, pch=16)
#' legend('topright', legend=c('retained', 'discarded'),
#' col=c('black', 'red'), pch=c(16, 1))
#' @export

geoThin <- function(
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

	# pairwise point distances
	dists <- geosphere::distm(xy, fun=distFunct)
	diag(dists) <- Inf

	# remove distances ("dists" will be changed as points are removed)
	masterDists <- dists

	# neighbors
	neighs <- colSums(dists <= minDist, na.rm=TRUE)

	if (any(neighs > 0)) {
	
		# while there is at least one point with neighbors too close
		if (verbose) {
			omnibus::say('Thinning points...')
			totalNeighs <- sum(neighs)
			prog <- utils::txtProgressBar(min=0, max=totalNeighs, width=32, style=3)
		}

		while (any(neighs > 0)) {

			# get points with greatest number of neighbors
			maxNeighs <- max(neighs)
			removeThis <- which(neighs == maxNeighs)

			# if more than one point has the maximum number of neighbors, break ties
			if (length(removeThis) > 1) {

				# removing point with closest neighbor
				distClosestNeigh <- apply(dists[ , removeThis], 2, min)
				minDistClosestNeigh <- min(distClosestNeigh)
				
				# catch cases where points are all at same coordinates, remove one randomly
				removeThis <- if (all(minDistClosestNeigh == 0)) {
					removeThis[sample(seq_along(removeThis), 1)]
				} else {
					removeThis[which(distClosestNeigh == minDistClosestNeigh)]
				}

				# if there is still a tie, break ties
				if (length(removeThis) > 1) {

					# break ties by removing point closest to centroid of all remaining points
					cent <- rgeos::gCentroid(xy[index])
					distToCent <- distFunct(cent, xy[removeThis], ...)
					minDistToCent <- min(distToCent)

					removeThis <- removeThis[which(distToCent == minDistToCent)]

					# if there is *still* a tie, break ties
					if (length(removeThis) > 1) {

						# break ties by removing point with least median distance to all other points (even if > minDist away)
						distToAllPoints <- apply(dists[ , removeThis], 2, median, na.rm=TRUE)
						minMedianDistToAllPoints <- min(distToAllPoints)
						removeThis <- removeThis[which(distToAllPoints == minMedianDistToAllPoints)]

						# if there is *still* a tie, break ties
						if (length(removeThis) > 1) {

							# break ties by randomly selecting a point
							removeThis <- removeThis[sample(seq_along(removeThis), 1)]

						}

					}

				}

			}

			# remove point index from those to keep
			index <- index[-which(index == removeThis)]
			
			dists[removeThis, ] <- rep(Inf, ncol(dists))
			dists[ , removeThis] <- rep(Inf, nrow(dists))

			# neighbors
			neighs <- colSums(dists <= minDist, na.rm=TRUE)

			if (verbose) utils::setTxtProgressBar(prog, totalNeighs - sum(neighs))

		}

		if (verbose) close(prog)

		# remove too-neighborly points
		x <- if (class(x) == 'data.frame' | class(x) == 'matrix') {
			x[index, ]
		} else if (class(x) == 'SpatialPointsDataFrame') {
			x[index, ]
		} else {
			x[index]
		}

	# all pairwise distances are > minimum distance
	} else {
	
		if (verbose) omnibus::say('No sites were within the minimum distance, so none were removed.')
		
	}

	x

}
