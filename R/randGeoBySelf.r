#' Randomizes geographic points while observing spatial autocorrelation
#'
#' This function randomizes the location of geographic points while retaining (more or less) the same distribution of pairwise distances between points (plus or minus a user-defined percent). The procedure is meant to generalized the "RTR" (rotate/translate/reflect) randomization procedure proposed by Nunes, L.A. and Pearson, R.G.  2017.  A null biogeographical test for assessing ecological niche evolution.  \emph{Journal of Biogeography} 44:1331-1343. The procedure in this function is actually adapted from the randomization procedure presented in Beale, C.M., J.J. Lennon, and A. Gimona.  2008.  Opening the climate envelope reveals no macroscale associations with climate in European birds. \emph{Proceedings of the National Academy of Sciences USA} 105:14908-14912.
#' @param x Matrix, data frame, SpatialPoints, or SpatialPointsDataFrame object. If this is a matrix or data frame, the first two columns must represent longitude and latitude (in that order). If \code{x} is a matrix or data frame, the coordinates are assumed to be unprojected (WGS84) (a coordinate reference system proj4 string or \code{CRS} object can be passed into the function using \code{...}). If \code{x} is a SpatialPoints or SpatialPointsDataFrame and not in WGS84 or NAD83, then coordinates are projected to WGS84 (with a warning).
#' @param rast Raster, RasterStack, or RasterBrick used to locate presences randomly. If this is a RasterStack or a RasterBrick then the first layer will be used (i.e., so cells with \code{NA} will not have points located within them).
#' @param bins Integer > 1, number of overlapping bins across which to calculate distribution of pairwise distances between points. The range covered by bins starts at 0 and and at the largest observed pairwise distance + 0.5 * bin width. Default value is 20.
#' @param tol Numeric >0, for any one bin, root-mean square deviation between observed pairwise distribution of distances and randomized distances required for the randomized distances to be considered statistically the same as observed distances.
#' @param distFunct Either a function or \code{NULL}. If \code{NULL} then \code{\link[geosphere]{distCosine}} is used to calculate distances.  Other "dist" functions (e.g., \code{\link[geosphere]{distGeo}}) can be used.  Alternatively, a custom function can be used so long as its first argument is a 2-column numeric matrix with one row for the x- and y-coordinates of a single point and its second argument is a two-column numeric matrix with one or more rows of other points.
#' @param verbose Logical, if \code{FALSE} (default) show no progress indicator. If \code{TRUE} then display occasional updates and graph.
#' @param ... Arguments to pass to \code{distCosine} or \code{\link[dismo]{randomPoints}}. Note that if \code{x} is a matrix or data frame a coordinate reference system may be passed using \code{crs = <proj4 string code} or \code{crs = <object of class CRS (see sp package)>}. Otherwise WGS84 is assumed.
#' @return Object of the same class as \code{x} but with coordinates randomized.
#' @seealso \code{\link[dismo]{randomPoints}}
#' @examples
#' # madagascar
#' data(mad)
#' par(layout(matrix(c(1, 2), nrow=1)))
#' plot(mad, main='Madagascar')
#' data(lemur)
#' points(lemur, pch=16)
#' rands <- randGeoBySelf(lemur, mad, verbose=TRUE)
#' par(fig=1, new=FALSE)
#' points(rand, col='red')
#' @export

randGeoBySelf <- function(
	x,
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
	out <- x

	# convert SpatialPointsDataFrame to SpatialPoints
	if (class(x) == 'SpatialPointsDataFrame') x <- sp::SpatialPoints(coordinates(x), proj4string=CRS(raster::projection(x)))

	# convert matrix/data frame to SpatialPoints
	if (class(x) %in% c('matrix', 'data.frame')) {

		x <- if (exists('crs', inherits=FALSE)) {
			sp::SpatialPoints(x[ , 1:2, drop=FALSE], sp::CRS(crs))
		} else {
			sp::SpatialPoints(x[ , 1:2, drop=FALSE], getCRS('wgs84', TRUE))
		}

	}

	# correct CRS
	if (raster::projection(x) != getCRS('wgs84') & raster::projection(x) != getCRS('nad83')) {

		warning('Coordinates are not in WGS84 or NAD83. Projecting them to WGS84.')
		x <- sp::spTransform(x, getCRS('wgs84', TRUE))

	}
	
	# remember CRS
	crs <- if ('crs' %in% omnibus::ellipseNames(list)) {
		ellipses$crs
	} else {
		raster::projection(x)
	}

	# check CRS of raster
	if (crs != raster::projection(rast)) {
		stop('Raster named in argument "rast" does not have same coordinate reference system as object named in "x" (or "crs").')
	}
	
	### calculate observed pairwise distances
	#########################################
	
	obsDists <- enmSdm::pointDist(x, ...)
	obsDists[upper.tri(obsDists, diag=TRUE)] <- NA
	obsDists <- c(obsDists)
	
	maxDist <- max(obsDists, na.rm=T)
	breaks <- c(0, maxDist * 1.1, bins)
	
	obsDistDistrib <- omnibus::histOverlap(obsDists, breaks=breaks)

	# randomize points: start by getting a large number... will cycle through these (faster than getting one-by-one)
	numRandPoints <- max(10000, length(x) * round(length(x) * bins / (10000 * tol)))
	randPoints <- enmSdm::sampleRast(rast, numRandPoints, prob=FALSE)
	randPoints <-  sp::SpatialPoints(randPoints, sp::CRS(crs))
	randSites <- randPoints[1:length(x)]
	randUsed <- length(x)
	
	randDists <- enmSdm::pointDist(randSites, ...)
	randDists[upper.tri(randDists, diag=TRUE)] <- NA
	randDists <- randDists
	
	randDistDistrib <- omnibus::histOverlap(randDists, breaks=breaks)

	# differences between observed and randomized distances
	delta <- sqrt(sum((randDistDistrib[ , 'proportion'] - obsDistDistrib[ , 'proportion'])^2))
	
	if (verbose) {

		par(mfrow=c(1, 1))
	
		plot(rast)
		points(x, pch=16)
		points(randSites, col='red')
		
		mids <- apply(obsDistDistrib[ , 1:2], 1, mean)
		plot(mids, obsDistDistrib[ , 'proportion'], pch=16, type='b', ylab='Proportion of Pairwise Distances', xlab='Distance Bin Midpoint')
		lines(mids, randDistDistrib[ , 'proportion'], col='red')
		legend('topright', inset=0.01, legend=c('Observed', 'Randomized'), col=c('black', 'red'), pch=c(16, NA), lwd=1)
		
	}
		
	tries <- accepts <- 0
		
	# iteratively randomized, check to see if this made distribution of randomized distances closer to observed distribution, if so keep
	while (delta > tol) {
	
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
		
		replaceIndex <- sample(seq_along(randSites), 1)
		candDist <- enmSdm::pointDist(candSite, randSites)
		candDist <- c(candDist)
		candRandDists <- randDists
		candRandDists[replaceIndex, ] <- candRandDists[ , replaceIndex] <- candDist
		candRandDists[upper.tri(candRandDists, diag=TRUE)] <- NA

		candRandDistDistrib <- omnibus::histOverlap(c(candRandDists), breaks=breaks)

		candDelta <- sqrt(sum((candDistDistrib[ , 'proportion'] - obsDistDistrib[ , 'proportion'])^2))
		
		# accept randomized point
		if (candDelta < delta | tries %% 10000 == 0) {

			if (verbose) {
				
				lines(mids, randDistDistrib[ , 'proportion'], col='gray80')
				lines(mids, candDistDistrib[ , 'proportion'], col='red')
				
			}
			
			coords <- sp::coordinates(randSites)
			coords[replaceIndex, ] <- sp::coordinates(candSite)
			randSites <- sp::SpatialPoints(coords, CRS(crs))
			
			randDistDistrib <- candDistDistrib
			randDists <- candRandDists
			accepts <- accepts + 1
			delta <- candDelta
				
			if (verbose) say('current tolerance: ', sprintf('%.6f', delta), ' | accepted: ', accepts, ' of ', tries, ' tries.')
			
		}
		
	}
	
	coords <- sp::coordinates(randSites)
	
	if (class(out) == 'SpatialPointsDataFrame') {
		out <- sp::SpatialPointsDataFrame(coords, data=as.data.frame(out), CRS(crs))
	} else if (class(out) == 'SpatialPoints') {
		out <- sp::SpatialPoints(coords, CRS(crs))
	} else {
		out[ , 1:2] <- coords
	}
	
	out	

}
