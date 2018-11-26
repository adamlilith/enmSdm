#' Randomizes geographic points while observing spatial autocorrelation
#'
#' This function randomizes the location of geographic points while retaining (more or less) the same distribution of pairwise distances between points (plus or minus a user-defined tolerance). The procedure is meant to generalize the "RTR" (rotate/translate/reflect) randomization procedure proposed by Nunes, L.A. and Pearson, R.G. 2017. A null biogeographical test for assessing ecological niche evolution. \emph{Journal of Biogeography} 44:1331-1343. The algorithm is a conceptual adaptation of the randomization procedure presented in Beale, C.M., J.J. Lennon, and A. Gimona. 2008. Opening the climate envelope reveals no macroscale associations with climate in European birds. \emph{Proceedings of the National Academy of Sciences USA} 105:14908-14912.
#' @param x Matrix, data frame, SpatialPoints, or SpatialPointsDataFrame object. If this is a matrix or data frame, the first two columns must represent longitude and latitude (in that order). If \code{x} is a matrix or data frame, the coordinates are assumed to be unprojected (WGS84) (a coordinate reference system proj4 string or \code{CRS} object can be passed into the function using \code{...}). If \code{x} is a SpatialPoints or SpatialPointsDataFrame and not in WGS84 or NAD83, then coordinates are projected to WGS84 (with a warning).
#' @param rast Raster, RasterStack, or RasterBrick demarcating the area in which randomized sites are to be placed. If a RasterStack or a RasterBrick is used then the first layer will be used (i.e., so cells with \code{NA} will not have points located within them).
#' @param tol Numeric >0, maximum root-mean-square distance allowed between the set of observed pairwise distances between points and the set of randomized pairwise distances between points. The algorithm will shuffle points until the calculated difference is <= this number. Units are the same as units used by the coordinate reference system of \code{x} (usually meters).
#' @param w Either \code{NULL} (default), a numeric matrix, or a function. This argument can be used to specify a set of weights to apply to the root-mean-square distance matrix. It must either be \code{NULL}, in which case all squared distances have equal weight, a symmetrical matrix with the same number of rows and columns as the length of \code{x} (diagonal elements are ignored), or a function. If \code{w} is a matrix then they are applied as \code{sqrt(sum(w * (d_obs - d_rand)^2) / sum(w))} where \code{d_obs} is the observed pairwise distance between points and \code{d_rand} the pairwise distance between the randomized points. If it is a function it is applied as \code{sqrt(mean(w(d_obs - d_rand)^2))}. Note that any weights other than \code{NULL} will also affect the desirable value of \code{tol}. This argument is especially useful for determining the scale of structure desired to be replicated in the randomized data. For example, using an increasing function will cause the randomization to recreate broad-scale patterns in the data (e.g., general shape and extent of \code{x}) but allow fine-scale patterns to be more random. A decreasing function will allow the broad pattern of \code{x} to shift (e.g., "bend") but retain fine-scale pattern.
#' @param distFunct Either a function or \code{NULL}. If \code{NULL} then \code{\link[geosphere]{distCosine}} is used to calculate distances. Other "dist" functions (e.g., \code{\link[geosphere]{distGeo}}) can be used for lower/increased accuracy and faster/slower performance. Alternatively, a custom function can be used so long as its input and output is comparable to any of these functions.
#' @param restrict Logical, if \code{TRUE} (default), then select random points from any non-\code{NA} cell in a restricted area. If \code{FALSE}, then select from any non-\code{NA} cell on the raster. If \code{FALSE}, then random sites are selected from across the entire raster.
#' @param keepData Logical, if \code{TRUE} then the original data in \code{x} (i.e., columns that do not represent coordinates) will be retained in the output but the coordinates will be shuffled. If \code{FALSE} (default) then the returned value will have just shuffled coordinates.
#' @param verbose Logical, if \code{FALSE} (default) show no progress indicator. If \code{TRUE} then display progress.
#' @param ... Arguments to pass to \code{distCosine} or \code{\link[enmSdm]{sampleRast}}. Note that if \code{x} is a matrix or data frame a coordinate reference system may be passed using \code{crs = <proj4 string code>} or \code{crs = <object of class CRS>} (see \pkg{sp} package). Otherwise the coordinates are assumed to be unprojected (WGS84).
#' @return Object of the same class as \code{x} but with coordinates randomized.
#' @details The argument \code{startBy} can be used to speed convergence to a solution. The two options, \code{anywhere} and \code{restricted}, start with patterns that are very over-dispersed or somewhat over-dispersed, respectively, but the \code{restricted} option will be slower to start. Under this option start up is initiated by placing a single random point on the landscape. Then, a buffer with a width equal to 1.05 * half the average maximum pairwise distance is drawn around the point. The raster in \code{rast} is cropped to the extent of this buffer and rasterized to form a circle. Note that if the centroid of the circle lies toward the edge of the raster or if the raster contains \code{NA} values then the circle may be cropped or may otherwise have "holes". Then, a number of random points equal to the number represented by \code{x} are placed within this area. The procedure then starts using these points. Note that the restricted start-up area is only used for the initial placement of points--subsequent randomizations use the entire study region.
#' @seealso \code{\link[enmSdm]{sampleRast}}
#' @examples
#' \donttest{
#' data(lemur)
#' data(mad0)
#' elev <- raster::getData('alt', country='MDG')
#' # thin points so none are <= 5000 m apart
#' # (to save time and reduce sampling bias)
#' lemur <- geoThin(lemur, minDist=7500)
#' Note: We are using a large "tol" to speed up the example.
#' set.seed(123)
#' rands <- randPointsRespectingSelf(lemur, elev, tol=15000, verbose=TRUE)
#' windows()
#' plot(elev, main='Madagascar')
#' plot(mad0, add=TRUE)
#' points(lemur, pch=16)
#' points(rands, col='red')
#' legend('bottomright', legend=c('Observed', 'Randomized'),
#'  col=c('black', 'red'), pch=c(16, 1))
#' }
#' @export
randPointsRespectingSelf <- function(
	x,
	rast,
	tol = NULL,
	w = NULL,
	distFunct = NULL,
	restrict = TRUE,
	keepData = FALSE,
	verbose = FALSE,
	...
) {

	ellipses <- list(...)

	if (is.null(distFunct)) distFunct <- geosphere::distCosine

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

	if (!is.null(w)) {
		if (class(w) != 'function') {
			if (dim(w)[1] != length(x) | dim(w)[2] != length(x)) {
				stop('Argument "w" must be a function or a matrix with number\nof rows and columns equal to the length of "x".')
			}
		}
	}
	
	# remember CRS
	crs <- if ('crs' %in% ellipses) {
		ellipses$crs
	} else {
		raster::projection(rast)
	}

	### calculate observed pairwise distances
	#########################################

	obsDists <- geosphere::distm(x, fun=distFunct)
	obsDists[upper.tri(obsDists, diag=TRUE)] <- NA
	obsDistsWithVals <- 2:nrow(obsDists)
	
	xIndex <- seq_along(x)
	
	if (is.null(tol)) {
	
		tol <- .calcTol(obsDists, TRUE)
		if (verbose) omnibus::say('Automatically calculated value for tol:  ', sprintf('%.2f', tol), '.')
		
	}

	### initiate random points
	##########################

	if (verbose) omnibus::say('Locating set candidate sites...')

	numRandPoints <- 10 * max(10000, round(length(x)^2 / max(1, log(tol))))

	# circle buffer width
	maxPairwiseDists <- apply(obsDists[obsDistsWithVals, ], 1, max, na.rm=TRUE)
	bufferWidth <- 0.55 * max(maxPairwiseDists)
	
	samples <- .sampleRastSubregion(
		rast = rast,
		n = numRandPoints,
		width = bufferWidth,
		restrict = restrict,
		circleCent = NULL,
		...
	)
	
	randPoints <- samples$randPoints
	circleCent <- samples$circleCent
	
	randSites <- randPoints[xIndex]
	randUsed <- length(x)
	
	randDists <- geosphere::distm(randSites, fun=distFunct)
	randDists[upper.tri(randDists, diag=TRUE)] <- NA

	delta <- statisfactory::rmsd(obsDists, randDists, w=w, na.rm=TRUE)
	
	### iteratively randomize
	#########################
	
	tries <- accepts <- sets <- 0

	# iteratively randomized, check to see if this made distribution of randomized distances closer to observed distribution, if so keep
	while (delta > tol) {

		tries <- tries + 1
		randUsed <- randUsed + 1

		# get new random set of coordinates
		if (randUsed > length(randPoints)) {

			if (verbose) omnibus::say('Locating another set of candidate sites (', sets + 1, ' so far', ')...')
			
			bufferWidth <- 1.05 * bufferWidth
			
			samples <- .sampleRastSubregion(
				rast = rast,
				n = numRandPoints,
				width = bufferWidth,
				restrict = TRUE,
				circleCent = circleCent,
				...
			)

			randPoints <- samples$randPoints
			randUsed <- 1
		
		}

		# replace selected random site with candidate random coordinate
		candSite <- randPoints[randUsed]

		replaceIndex <- sample(xIndex, 1)
		candDist <- geosphere::distm(candSite, randSites, fun=distFunct)
		
		candRandDists <- randDists
		candRandDists[replaceIndex, ] <- candRandDists[ , replaceIndex] <- candDist
		candRandDists[upper.tri(candRandDists, diag=TRUE)] <- NA

		candDelta <- statisfactory::rmsd(obsDists, candRandDists, w=w, na.rm=TRUE)
		
		# accept randomized point
		if (candDelta < delta) {

			if (verbose) {
				omnibus::say('(actual - desired) tolerance: ', sprintf('%.2f', candDelta - tol), ' with ', accepts, ' accepts of ', tries, ' tries')
			}

			randSites@coords[replaceIndex, ] <- sp::coordinates(candSite)

			randDists <- candRandDists
			accepts <- accepts + 1
			delta <- candDelta

		}

	}

	out <- .randPointsReturn(x=x, randSites=randSites, crs=crs, keepData=keepData, out=out) 

	attr(out, 'keepData') <- keepData
	attr(out, 'weights') <- w
	attr(out, 'distFunct') <- distFunct
	attr(out, 'tries') <- tries
	attr(out, 'accepts') <- accepts
	attr(out, 'restrict') <- restrict
	
	attr(out, 'tol') <- tol
	attr(out, 'finalTol') <- delta

	out

}
