#' Randomizes the location of two sets of geographic points while respecting spatial autocorrelation
#'
#' This function randomizes the location of two sets of geographic points with respect to one another retaining (more or less) the same distribution of pairwise distances between points with and between sets (plus or minus a user-defined tolerance).
#' @param x1 Matrix, data frame, SpatialPoints, or SpatialPointsDataFrame object. If this is a matrix or data frame, the first two columns must represent longitude and latitude (in that order). If \code{x} is a matrix or data frame, the coordinates are assumed to be unprojected (WGS84) (a coordinate reference system proj4 string or \code{CRS} object can be passed into the function using \code{...}). If \code{x} is a SpatialPoints or SpatialPointsDataFrame and not in WGS84 or NAD83, then coordinates are projected to WGS84 (with a warning).
#' @param x2 As \code{x1}.
#' @param rast Raster, RasterStack, or RasterBrick used to locate presences randomly. If this is a RasterStack or a RasterBrick then the first layer will be used (i.e., so cells with \code{NA} will not have points located within them).
#' @param bins Integer > 1, number of overlapping bins across which to calculate distribution of pairwise distances between points. The range covered by bins starts at 0 and and at the largest observed pairwise distance + 0.5 * bin width. Default value is 20.
#' @param tol Numeric >0, for any one bin, root-mean square deviation between observed pairwise distribution of distances and randomized distances required for the randomized distances to be considered statistically the same as observed distances.
#' @param distFunct Either a function or \code{NULL}. If \code{NULL} then \code{\link[geosphere]{distCosine}} is used to calculate distances.  Other "dist" functions (e.g., \code{\link[geosphere]{distGeo}}) can be used.  Alternatively, a custom function can be used so long as its first argument is a 2-column numeric matrix with one row for the x- and y-coordinates of a single point and its second argument is a two-column numeric matrix with one or more rows of other points.
#' @param verbose Logical, if \code{FALSE} (default) show no progress indicator. If \code{TRUE} then display occasional updates and graph.
#' @param ... Arguments to pass to \code{distCosine} or \code{\link[dismo]{randomPoints}}. Note that if \code{x} is a matrix or data frame a coordinate reference system may be passed using \code{crs = <proj4 string code} or \code{crs = <object of class CRS (see sp package)>}. Otherwise WGS84 is assumed.
#' @return Object of the same class as \code{x} but with coordinates randomized.
#' @seealso \code{\link[dismo]{randomPoints}}, \code{\link[enmSdm]{randGeoBySelf}}
#' @examples
#' # madagascar
#' library(dismo)
#' madElev <- getData('alt', country='MDG')
#' par(layout(matrix(c(1, 2), nrow=1)))
#' plot(madElev, main='Madagascar')
#' data(lemur)
#' points(lemur, pch=16)
#' rands <- randGeoBySelf(lemur, mad, verbose=TRUE)
#' par(fig=1, new=FALSE)
#' points(rand, col='red')
#' @export

randGeoBySelfAndOther <- function(
	x1,
	x2,
	rast,
	bins=20,
	tol=0.001,
	distFunct=NULL,
	verbose=FALSE,
	...
) {

	if (is.null(distFunct)) distFunct <- geosphere::distCosine

	ellipses <- list(...)

	# save copy of original
	out1 <- x1
	out2 <- x2

	# convert SpatialPointsDataFrame to SpatialPoints
	if (class(x1) == 'SpatialPointsDataFrame') x1 <- sp::SpatialPoints(coordinates(x1), proj4string=CRS(raster::projection(x1)))
	if (class(x2) == 'SpatialPointsDataFrame') x2 <- sp::SpatialPoints(coordinates(x2), proj4string=CRS(raster::projection(x2)))

	# convert matrix/data frame to SpatialPoints
	if (class(x1) %in% c('matrix', 'data.frame')) {

		x1 <- if (exists('crs', inherits=FALSE)) {
			sp::SpatialPoints(x1[ , 1:2, drop=FALSE], sp::CRS(crs))
		} else {
			sp::SpatialPoints(x1[ , 1:2, drop=FALSE], getCRS('wgs84', TRUE))
		}

	}

	if (class(x2) %in% c('matrix', 'data.frame')) {

		x2 <- if (exists('crs', inherits=FALSE)) {
			sp::SpatialPoints(x2[ , 1:2, drop=FALSE], sp::CRS(crs))
		} else {
			sp::SpatialPoints(x2[ , 1:2, drop=FALSE], getCRS('wgs84', TRUE))
		}

	}

	# correct CRS
	if (!(raster::projection(x1) == getCRS('wgs84') | raster::projection(x1) == getCRS('nad83'))) {

		warning('Coordinates of x1 are not in WGS84 or NAD83 coordinate reference system. Projecting them to WGS84.')
		x1 <- sp::spTransform(x1, getCRS('wgs84', TRUE))

	}
	
	if (!(raster::projection(x2) != getCRS('wgs84') | raster::projection(x2) != getCRS('nad83'))) {

		warning('Coordinates of x2 are not in WGS84 or NAD83 coordinate reference system. Projecting them to WGS84.')
		x2 <- sp::spTransform(x2, getCRS('wgs84', TRUE))

	}

	if (raster::projection(x1) != raster::projection(x2)) stop('Coordinate reference systems for x1 and x2 are not the same.')
	
	# remember CRS
	crs <- if ('crs' %in% omnibus::ellipseNames(list)) {
		ellipses$crs
	} else {
		raster::projection(x1)
	}

	# check CRS of raster
	if (crs != raster::projection(rast)) {
		stop('Raster named in argument "rast" does not have same coordinate reference system as objects named in "x1" or "x2".')
	}
	
	### calculate observed pairwise distances within and between point sets
	#######################################################################
	
	obsSelfDists1 <- enmSdm::pointDist(x1, ...)
	obsSelfDists1[upper.tri(obsSelfDists1, diag=TRUE)] <- NA
	obsSelfDists1 <- c(obsSelfDists1)
	
	obsSelfDists2 <- enmSdm::pointDist(x2, ...)
	obsSelfDists2[upper.tri(obsSelfDists2, diag=TRUE)] <- NA
	obsSelfDists2 <- c(obsSelfDists2)

	obsOtherDists <- enmSdm::pointDist(x1, x2, ...)
	for (i in 2:(nrow(obsOtherDists) - 1)) obsOtherDists[i, (i + 1):ncol(obsOtherDists)] <- NA
	
	maxSelfDist1 <- max(obsSelfDists1, na.rm=TRUE)
	breaksSelf1 <- c(0, maxSelfDist1 * 1.05, bins)
	
	maxSelfDist2 <- max(obsSelfDists2, na.rm=TRUE)
	breaksSelf2 <- c(0, maxSelfDist2 * 1.05, bins)
	
	maxOtherDist <- max(obsOtherDists, na.rm=TRUE)
	breaksOther <- c(0, maxOtherDist * 1.05, bins)
	
	obsSelfDistDistrib1 <- omnibus::histOverlap(obsSelfDists1, breaks=breaksSelf1)
	obsSelfDistDistrib2 <- omnibus::histOverlap(obsSelfDists2, breaks=breaksSelf2)
	obsOtherDistDistrib <- omnibus::histOverlap(obsOtherDists, breaks=breaksOther)

	# randomize points: start by getting a large number...
	# will cycle through these (faster than getting one-by-one)
	x1size <- length(x1)
	x2size <- length(x2)
	x12size <- x1size + x2size
	
	numRandPoints <- max(10000, x12size * round(x12size * bins / (10000 * tol)))
	
	randPoints1 <- enmSdm::sampleRast(rast, numRandPoints, prob=FALSE)
	randPoints2 <- enmSdm::sampleRast(rast, numRandPoints, prob=FALSE)
	
	randPoints1 <-  sp::SpatialPoints(randPoints1, sp::CRS(crs))
	randPoints2 <-  sp::SpatialPoints(randPoints2, sp::CRS(crs))
	
	# allocate randomized points to each set
	randSites1 <- randPoints[1:x1size]
	randSites2 <- randPoints[(x1size + 1):x12size]
	
	randUsed <- x12size
	
	# calculate intra- and inter-set pairwise distances
	randSelfDists1 <- enmSdm::pointDist(randSites1, distFunct=distFunct, ...)
	randSelfDists2 <- enmSdm::pointDist(randSites2, distFunct=distFunct, ...)
	randSelfDists1[upper.tri(randDists1, diag=TRUE)] <- NA
	randSelfDists2[upper.tri(randDists2, diag=TRUE)] <- NA
	
	randOtherDists <- enmSdm::pointDist(randSites1, randSites2, distFunct=distFunct, ...)
	for (i in 2:(nrow(obsOtherDists) - 1)) {
		obsOtherDists[i, (i + 1):ncol(obsOtherDists)] <- NA
	}
	
	randSelfDists1 <- c(randSelfDists1)
	randSelfDists2 <- c(randSelfDists2)
	randOtherDists <- c(randOtherDists)
	
	randSelfDistDistrib1 <- omnibus::histOverlap(randSelfDists1, breaks=breaksSelf1)
	randSelfDistDistrib2 <- omnibus::histOverlap(randSelfDists2, breaks=breaksSelf2)
	randOtherDistDistrib <- omnibus::histOverlap(randOtherDists, breaks=breaksOther)

	# differences between observed and randomized distances
	deltaSelf1 <- sqrt(sum((randSelfDistDistrib1[ , 'proportion'] - obsSelfDistDistrib1[ , 'proportion'])^2))
	deltaSelf2 <- sqrt(sum((randSelfDistDistrib2[ , 'proportion'] - obsSelfDistDistrib2[ , 'proportion'])^2))
	deltaOther <- sqrt(sum((randOtherDistDistrib[ , 'proportion'] - obsOtherDistDistrib[ , 'proportion'])^2))

	# plot
	if (verbose) {

		par(mfrow=c(1, 4))
	
		plot(rast)
		points(x1, pch=1)
		points(x2, pch=2)
		points(randSites1, col='red', pch=1)
		points(randSites2, col='blue', pch=2)
		
		# x1
		midsSelf1 <- apply(obsSelfDistDistrib1[ , 1:2], 1, mean)
		plot(midsSelf1, obsSelfDistDistrib1[ , 'proportion'], pch=16, type='b', ylab='Proportion of Pairwise Distances (x1)', xlab='Distance Bin Midpoint')
		lines(midsSelf1, randSelfDistDistrib1[ , 'proportion'], col='red')
		legend('topright', inset=0.01, legend=c('Observed', 'Randomized'), col=c('black', 'red'), pch=c(16, NA), lwd=1)
		
		# x2
		midsSelf2 <- apply(obsSelfDistDistrib2[ , 1:2], 1, mean)
		plot(midsSelf2, obsSelfDistDistrib2[ , 'proportion'], pch=16, type='b', ylab='Proportion of Pairwise Distances (x2)', xlab='Distance Bin Midpoint')
		lines(midsSelf2, randSelfDistDistrib2[ , 'proportion'], col='red')
		legend('topright', inset=0.01, legend=c('Observed', 'Randomized'), col=c('black', 'red'), pch=c(16, NA), lwd=1)
		
		# vs other
		midsOther <- apply(obsOtherDistDistrib[ , 1:2], 1, mean)
		plot(midsOther, obsOtherDistDistrib[ , 'proportion'], pch=16, type='b', ylab='Proportion of Pairwise Distances (x2 vs x2)', xlab='Distance Bin Midpoint')
		lines(midsOther, randOtherDistDistrib[ , 'proportion'], col='red')
		legend('topright', inset=0.01, legend=c('Observed', 'Randomized'), col=c('black', 'red'), pch=c(16, NA), lwd=1)
		
	}
		
	tries <- accepts <- 0
		
	# iteratively randomized, check to see if this made distribution of randomized distances closer to observed distribution, if so keep
	while (deltaSelf1 > tol && deltaSelf2 > tol && deltaOther > tol) {
	
		tries <- tries + 1
		randUsed <- randUsed + 1
	
		# get new random set of coordinates
		if (randUsed > length(randPoints)) {
			randPoints <- enmSdm::sampleRast(rast, numRandPoints, prob=FALSE)
			randPoints <- sp::SpatialPoints(randPoints, sp::CRS(crs))
			randUsed <- 1
		}

		# replace selected random site with candidate random coordinate
		candSite <- randPoints[randUsed]

		# randomly choose among sets
		whichSet <- runif(1)
		
		chosenSet <- if (whichSet <= x1size / x12size) {
			'x1'
		} else {
			'x2'
		}
		
		### try replacing one point in set 1
		if (chosenSet == 'x1') {

			replaceIndex <- sample(seq_along(randSites1), 1)
			candDistVsSelf <- enmSdm::pointDist(candSite, randSites1)
			
			# distribution of distances vs self
			candRandDists1 <- randDists1
			candRandDists1[replaceIndex, ] <- candRandDists1[ , replaceIndex] <- candDistVsSelf
			candRandDists1[upper.tri(candRandDists1, diag=TRUE)] <- NA
			candRandDistDistrib1 <- omnibus::histOverlap(c(candRandDists1), breaks=breaksSelf1)
			
			# distribution of distances of other
			candRandDistDistrib2 <- if (tries == 1) {
				omnibus::histOverlap(c(randSelfDists2), breaks=breaksSelf2)
			} else {
				randDistDistrib2
			}

			# distribution of distances vs other
			candDistVsOther <- randDistVsOther
			thisCandDistVsOther <- enmSdm::pointDist(candSite, randSites2, distFunct=distFunct, ...)
			candDistVsOther[replaceIndex, ] <- thisCandDistVsOther

			# differences between observed and expected
			candDeltaSelf1 <- sqrt(sum((candSelfDistDistrib1[ , 'proportion'] - obsSelfDistDistrib1[ , 'proportion'])^2))
			candDeltaSelf2 <- if (tries == 1) {
				sqrt(sum((randSelfDistDistrib2[ , 'proportion'] - obsSelfDistDistrib2[ , 'proportion'])^2))
			} else {
				randDeltaSelf2
			}

		### try replacing one point in set 2		
		} else {

			replaceIndex <- sample(seq_along(randSites2), 1)
			candDistVsSelf <- enmSdm::pointDist(candSite, randSites2)
			
			# distribution of distances vs self
			candRandDists2 <- randDists2
			candRandDists2[replaceIndex, ] <- candRandDists2[ , replaceIndex] <- candDistVsSelf
			candRandDists2[upper.tri(candRandDists2, diag=TRUE)] <- NA
			candRandDistDistrib2 <- omnibus::histOverlap(c(candRandDists2), breaks=breaksSelf2)
			
			# distribution of distances of other
			candRandDistDistrib1 <- if (tries == 1) {
				omnibus::histOverlap(c(randSelfDists1), breaks=breaksSelf1)
			} else {
				randDistDistrib1
			}

			# distribution of distances vs other
			candDistVsOther <- randDistVsOther
			thisCandDistVsOther <- enmSdm::pointDist(candSite, randSites1, distFunct=distFunct, ...)
			candDistVsOther[ , replaceIndex] <- thisCandDistVsOther

			# differences between observed and expected
			candDeltaSelf1 <- if (tries == 1) {
				sqrt(sum((candSelfDistDistrib1[ , 'proportion'] - obsSelfDistDistrib1[ , 'proportion'])^2))
			} else {
				randDeltaSelf1
			}
			candDeltaSelf2 <- sqrt(sum((randSelfDistDistrib2[ , 'proportion'] - obsSelfDistDistrib2[ , 'proportion'])^2))

		}
		
		# differences between observed and randomized distances
		candDeltaOther <- sqrt(sum((candOtherDistDistrib[ , 'proportion'] - obsOtherDistDistrib[ , 'proportion'])^2))
		
		# accept randomized point
		if (candDeltaSelf1 < deltaSelf1 & candDeltaSelf2 < deltaSelf2 & candDeltaOther < deltaOther | tries %% 10000 == 0) {

			if (verbose) {
				
				if (chosenSet == 'x1') {
					lines(mids1, candDistDistrib1[ , 'proportion'], col='gray80')
					lines(mids1, candDistDistrib1[ , 'proportion'], col='red')
				} else {
					lines(mids2, candDistDistrib2[ , 'proportion'], col='gray80')
					lines(mids2, candDistDistrib2[ , 'proportion'], col='red')
				}

			}
			
			if (chosenSet == 'x1') {
				
				coords <- sp::coordinates(randSites1)
				coords[replaceIndex, ] <- sp::coordinates(candSite)
				randSites1 <- sp::SpatialPoints(coords, CRS(crs))
				
				randDists1 <- candRandDists1
				candRandDistDistrib1 <- randDistDistrib1
				
			} else {
				
				coords <- sp::coordinates(randSites2)
				coords[replaceIndex, ] <- sp::coordinates(candSite)
				randSites2 <- sp::SpatialPoints(coords, CRS(crs))
				
				randDists2 <- candRandDists2
				candRandDistDistrib2 <- randDistDistrib2
				
			}
			
			randDistVsOther <- candDistVsOther
			
			accepts <- accepts + 1
			deltaSelf1 <- candDelta1
			deltaSelf2 <- candDelta2
			deltaOther <- candDeltaOther
				
			if (verbose) say('current tolerances: ', sprintf('%.6f', deltaSelf1), ' ', sprintf('%.6f', deltaSelf2), sprintf('%.6f', deltaOther), ' | accepted: ', accepts, ' of ', tries, ' tries.')
			
		}
		
	}
	
	coords1 <- sp::coordinates(randSites1)
	coords2 <- sp::coordinates(randSites2)
	
	if (class(out1) == 'SpatialPointsDataFrame') {
		out1 <- sp::SpatialPointsDataFrame(coords1, data=as.data.frame(out1), CRS(crs))
	} else if (class(out1) == 'SpatialPoints') {
		out1 <- sp::SpatialPoints(coords1, CRS(crs))
	} else {
		out1[ , 1:2] <- coords1
	}
	
	if (class(out2) == 'SpatialPointsDataFrame') {
		out2 <- sp::SpatialPointsDataFrame(coords2, data=as.data.frame(out2), CRS(crs))
	} else if (class(out2) == 'SpatialPoints') {
		out2 <- sp::SpatialPoints(coords2, CRS(crs))
	} else {
		out2[ , 1:2] <- coords2
	}
	
	out	<- list()
	out[[1]] <- out1
	out[[2]] <- out2
	out

}
