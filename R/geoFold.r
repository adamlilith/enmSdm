#' Assign geographically-distinct k-folds
#'
#' This function assigns geographically-divided k-folds ("g-folds") using partitioning around mediods (PAM) algorithm. The user can specify the number of folds to create, and optionally, the minimum size of any fold plus the minimum number of sites NOT in any fold (good for ensuring each fold has enough sites for testing and training).
#' @param x A data frame, matrix, \code{SpatialPoints}, or \code{SpatialPointsDataFrame} object. If a data frame or matrix then the coordinate reference system is assumed to be unprojected (WGS84).
#' @param k Positive integer. Number of k-folds to create.
#' @param minIn Positive integer or \code{NULL}. Minimum number of sites required to be in a fold. If left \code{NULL} (default), it is possible to have just one site in a fold.
#' @param minOut Positive integer or \code{NULL}. Minimum number of sites required to be outside of a fold. (i.e., if there are 5 folds, then for fold #1 this is the number of sites in folds 2 through 5). Leave as NULL to ignore.
#' @param longLat Two-element character list \emph{or} two-element integer list. If \code{x} is a data frame then this should be a character list specifiying the names of the fields in \code{x} \emph{or} a two-element list of integers that correspond to longitude and latitude (in that order). For example, \code{c('long', 'lat')} or \code{c(1, 2)}. If \code{x} is a matrix then this is a two-element list indicating the column numbers in \code{x} that represent longitude and latitude. For example, \code{c(1, 2)}. If \code{x} is a \code{SpatialPoints} object then this is ignored.
#' @param distFunct Either a function or \code{NULL}. If \code{NULL} then \code{distCosine()} in the \code{geoshere} package is used to calculate distances.  More accurate distances can be obtained by using other functions (see [geosphere::distCosine()] and related "\code{dist}" functions). Alternatively, a custom function can be used so long as its first argument is a 2-column numeric matrix with one row for the x- and y-coordinates of a single point and its second argument is a two-column numeric matrix with one or more rows of other points.
#' @param swaps Positive integer. Sometimes the routine generates folds that aren't minimally compact; i.e., points from some folds are spatially inside other folds. To correct this a random swap procedure is performed at the end in which pairs of points from different folds are swapped assignment. If this decreases the mean distance to the (new) centroid of each fold then the swap is kept. Otherwise it is not. This procedure is performed \code{swap} number of times.  Default (set by \code{swap=NULL}) is \code{max(100, n^2)} times where \code{n} is total number of points.
#' @param ... Arguments to pass to \code{distFunct}.
#' @return An integer vector, one element for for of x, with values 1 through k indicating which fold a site is located in.
#' @seealso \code{\link[geosphere]{distCosine}}, \code{\link[cluster]{pam}}
#' @examples
#' # Make three groups, one with one point and two with 20 points apiece.
#' # Naturally these should group into 3 groups with 1, 20, and 20 point apiece.
#' # By setting minIn and minOut to non-NULL values, we can increase/decrease
#' # the size of the groups.
#' # define plot function
#' pointPlot <- function(x, folds, ...) {
#'    plot(x, pch=16, cex=2, col='white', ...)
#'    for (i in sort(unique(folds))) points(x[folds==i, ], bg=i + 1, pch=20 + i, cex=2)
#'    legend('bottomright', legend=paste('fold', sort(unique(folds))),
#'    pt.bg=sort(unique(folds)) + 1, pch=20 + sort(unique(folds)), cex=1.4)
#' }
#'
#' set.seed(17)
#' group1 <- data.frame(x=c(-90, -90), y=c(40, 41))
#' group2 <- data.frame(x=rep(-80, 20), y=rep(37, 20))
#' group3 <- data.frame(x=rep(-100, 20), y=rep(37, 20))
#'
#' group2 <- group2 + cbind(rnorm(20), rnorm(20))
#' group3 <- group3 + cbind(rnorm(20), rnorm(20))
#'
#' sites <- rbind(group1, group2, group3)
#'
#' # simple g-folds
#' folds <- geoFold(sites, k=3)
#' pointPlot(sites, folds, main='Simple G-folds')
#'
#' # g-folds with minimum number of sites per fold
#' folds <- geoFold(sites, k=3, minIn=5)
#' pointPlot(sites, folds, main='G-folds with >=10 sites in each')
#'
#' # g-folds with minimum number of outside fold
#' folds <- geoFold(sites, k=3, minIn=10, minOut=10)
#' pointPlot(sites, folds, main='G-folds with >=10\nsites in/outside each')
#'
#' # g-folds with minimum number inside and of outside fold
#' folds <- geoFold(sites, k=3, minIn=14, minOut=20)
#' pointPlot(sites, folds, main='G-folds >=14 in\nand >=20 outside')
#' @export

geoFold <- function(
	x,
	k,
	minIn = NULL,
	minOut = NULL,
	longLat = NULL,
	distFunct = NULL,
	swaps = NULL,
	...
) {

	############################
	## functions and packages ##
	############################

	if (is.null(distFunct)) distFunct <- geosphere::distGeo

	# count number of sites outside each fold
	numOut <- function(folds) {

		out <- rep(NA, k)
		for (i in 1:k) out[i] <- sum(folds!=i, na.rm=T)
		out

	}

	# count number of sites outside each fold
	numIn <- function(folds) {

		inFold <- rep(NA, k)
		for (i in 1:k) inFold[i] <- sum(folds==i, na.rm=T)
		inFold

	}

	# get coordinates
	xy <- xToCoords(x, longLat)

	# catch impossible arguments... numIn and numOut too large
	if (!is.null(minIn) & !is.null(minOut)) if (k * minIn > length(xy)) stop('minOut must be <= k * minIn.')
	if (!is.null(minIn) & !is.null(minOut)) if (minOut > (k - 1) * minIn) stop('minOut must be <= (k - 1) * minIn.')
	if (!is.null(minIn) & !is.null(minOut)) if (length(xy) < minIn + minOut) stop('Total number of sites must be >= minIn + minOut.')

	if (is.null(minIn)) minIn <- 1
	if (is.null(minOut)) minOut <- 1

	##########
	## MAIN ##
	##########

	### get initial folds and test size
	###################################

	# generate folds
	dists <- geosphere::distm(xy, fun=distFunct)
	diag(dists) <- NA
	dists <- stats::as.dist(dists)
	folds <- as.integer(cluster::pam(x=dists, k=k, diss=TRUE, do.swap=T, cluster.only=T))

	shuffled <- FALSE # flag saying if points were shuffled between folds

	### shift sites to ensure minimum number in each fold
	#####################################################

	# general approach:
	# 1	find overall centroid C
	# 2 for each site calculate distance to C
	# 3 for each fold find site closest to C, call these c
	# 4 choose fold with c that is farthest from C
	# 5 calculate distance from each point in other folds to centroid of focal fold
	# 6 steal the point closest to the focal point's centroid
	# 7 repeat from 5 until focal fold has sufficient sites
	# 8 calculate new PAM using just points that do not belong to any former focal fold
	# 9 repeat from 1, excluding points that have been in focal folds

	if (!is.null(minIn)) {

		# counter for folds left to be handled
		kToDo <- 1:k

		while (any(numIn(folds) < minIn) | length(unique(folds)) < k) {

			shuffled <- TRUE

			origfolds <- folds

			# generate new folds for any sites that are not part of folds that have already been processed...
			# tempFolds contains new fold names PLUS NAs for folds that have been processed
			tempFolds <- folds
			interPointDists <- geosphere::distm(xy, fun=distFunct)
			diag(interPointDists) <- NA
			interPointDists <- interPointDists[which(folds %in% kToDo), which(folds %in% kToDo)]
			interPointDists <- as.dist(interPointDists)
			tempFolds[!(tempFolds %in% kToDo)] <- NA
			tempFolds[!is.na(tempFolds)] <- as.integer(cluster::pam(x=interPointDists, k=length(kToDo), diss=TRUE, do.swap=TRUE, cluster.only=TRUE))

			# increase fold counters by those that have been done
			for (i in (1:k)[-kToDo]) tempFolds[which(!is.na(tempFolds) & tempFolds >= i)] <- tempFolds[which(!is.na(tempFolds) & tempFolds >= i)] + 1

			# re-assign folds using names of folds already done plus new fold names from pam() above
			folds[which(folds %in% kToDo)] <- tempFolds[!is.na(tempFolds)]

			# find overall centroid of folds that can be manipulated
			centroid <- rgeos::gCentroid(xy[folds %in% kToDo])

			# distances to centroid for all points
			# dists <- distFunct(xy, centroid, ...)
			dists <- distFunct(xy, centroid)
			if (length(kToDo) < k) dists[which(!(folds %in% kToDo))] <- Inf

			# find minimum distance to centroid across folds
			closest <- rep(-Inf, k)
			for (i in kToDo) closest[i] <- min(dists[folds==i], na.rm=T)
			closest[numIn(folds) >= minIn] <- -Inf

			# get fold with farthest minimum distance that also needs more sites added to it
			focal <- which.max(closest)

			# for each point to steal
			for (countPoints in 1:(minIn - sum(folds==focal, na.rm=TRUE))) {

				# get distance between focal fold's centroid and all candidate points
				# dists <- distFunct(rgeos::gCentroid(xy[folds==focal]), xy, ...)
				dists <- distFunct(rgeos::gCentroid(xy[folds==focal]), xy)
				dists[folds==focal] <- Inf
				if (length(kToDo) < k) dists[!(folds %in% kToDo)] <- Inf

				# get index of point to steal and steal it
				stealFrom <- which.min(dists)
				folds[stealFrom] <- focal

			}

			kToDo <- kToDo[-which(kToDo==focal)]

		} # while there is at least one fold too small

	} # if minIn specified

	### shift sites to ensure minimum number outside of each fold
	#############################################################

	# general approach:
	# 1 find largest fold
	# 2 calculate distances from each point in largest to centroid of all other points
	# 3 find shortest distance between largests' points and other's centroids
	# 4 donate a point to another folds based on minimum distance to others' centroids
	# 5 repeat from 2 until focal fold has sufficent points outside it
	# 6 repeat from 1 but do not use folds that have already been focii

	# NB for shorthand, calling folds not part of a fold "outfolds" (e.g., fold #2's outfolds are 1, 3, 4, and 5 if k=5)

		if (!is.null(minOut)) {

			# counter for folds left to be handeled
			kToDo <- 1:k

			while (any(numOut(folds) < minOut) & length(kToDo) > 0) {

				shuffled <- TRUE

				# identify fold to donate from ... break ties by closeness to total centroid
				if (sum(numOut(folds)[kToDo]==min(numOut(folds)[kToDo])) == 1) {

					focal <- as.integer(names(which.max(numIn(folds)[kToDo])))

				} else {

					centroid <- rgeos::gCentroid(xy)
					dists <- rep(Inf, k)
					for (i in kToDo) dists[i] <- distFunct(rgeos::gCentroid(xy[folds==i]), centroid, ...)
					dists[numOut(folds)!=min(numOut(folds))] <- Inf

					focal <- which.min(dists)

				}

				# for each point to donate
				for (countPoints in 1:(minOut - numOut(folds)[focal])) {

					# get distance between focal fold's points and centroids of other eligible folds
					cents <- data.frame()
					for (i in 1:k) cents <- rbind(cents, sp::coordinates(rgeos::gCentroid(xy[folds==i])))

					dists <- matrix(rep(Inf, length(xy) * k), nrow=length(xy))
					for (i in which(folds==focal)) dists[i, ] <- distFunct(xy[i], cents, ...)
					if (length(kToDo) < k) dists[ , (1:k)[-kToDo]] <- Inf
					dists[ , focal] <- Inf

					minDistEachPoint <- apply(dists, 1, min)

					# point to donate
					donateFromIndex <- which.min(minDistEachPoint)

					# donate
					foldToDonateTo <- which.min(dists[donateFromIndex, ])
					folds[donateFromIndex] <- foldToDonateTo

				}

				kToDo <- kToDo[-which(kToDo==focal)]

			} # while any outfolds too small

		} # if minOut specified

	### randomly swap sites between folds to reduce spatial overlap between them
	############################################################################

	# general procedure:
	# 1 randomly choose any two points in different folds
	# 2 calculate mean distance from each fold's points to its centorid
	# 3 temporarily swap fold labels
	# 4 calculate mean distance from each fold's points to its new centorid with swapped labels
	# 5 if mean distance decreases for both then keep swapped labels, if not discard

	if (shuffled) {

		if (is.null(swaps)) swaps <- max(100, length(xy)^2)

		origfolds <- folds
		for (i in 1:swaps) {

			# get folds
			fold1 <- sample(1:k, 1)
			fold2 <- sample((1:k)[-fold1], 1)

			# swap labels
			tempFolds <- folds
			index1 <- sample(which(folds==fold1), 1)
			index2 <- sample(which(folds==fold2), 1)

			tempFolds[index1] <- fold2
			tempFolds[index2] <- fold1

			# unadulterated fold centroids
			centBefore1 <- rgeos::gCentroid(xy[folds==fold1])
			centBefore2 <- rgeos::gCentroid(xy[folds==fold2])

			# adulturated centroids
			centAfter1 <- rgeos::gCentroid(xy[tempFolds==fold1])
			centAfter2 <- rgeos::gCentroid(xy[tempFolds==fold2])

			# compare mean distances
			meanDistBefore1 <- mean(distFunct(centBefore1, xy[folds==fold1], ...))
			meanDistBefore2 <- mean(distFunct(centBefore2, xy[folds==fold2], ...))

			meanDistAfter1 <- mean(distFunct(centAfter1, xy[tempFolds==fold1], ...))
			meanDistAfter2 <- mean(distFunct(centAfter2, xy[tempFolds==fold2], ...))

			if (meanDistBefore1 > meanDistAfter1 & meanDistBefore2 > meanDistAfter2) folds <- tempFolds

		} # next swap

	}

	#####################
	## post-processing ##
	#####################

	folds

}
