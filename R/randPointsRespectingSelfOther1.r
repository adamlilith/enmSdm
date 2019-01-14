#' Randomizes the location of two sets of geographic points while respecting spatial autocorrelation
#'
#' THIS FUNCTION PROBABLY WONT WORK AS-IS--IT NEEDS A NEW RANDOMIZATION PROCEDURE TO ACHIEVE CONVERGENCE. This function randomizes the location of two sets of geographic points with respect to one another retaining (more or less) the same distribution of pairwise distances between points with and between sets (plus or minus a user-defined tolerance).
#' @param x1 Matrix, data frame, SpatialPoints, or SpatialPointsDataFrame object representing geographic points to be shuffled. If this is a matrix or data frame, the first two columns must represent longitude and latitude (in that order). If \code{x} is a matrix or data frame, the coordinates are assumed to be unprojected (WGS84) (a coordinate reference system proj4 string or \code{CRS} object can be passed into the function using \code{...}). If \code{x} is a SpatialPoints or SpatialPointsDataFrame and not in WGS84 or NAD83, then coordinates are projected to WGS84 (with a warning).
#' @param x2 As \code{x1}. Note that only points in \code{x1}. Points in \code{x2} are kept fixed on the landscape.
#' @param rast Raster, RasterStack, or RasterBrick used to locate presences randomly. If this is a RasterStack or a RasterBrick then the first layer will be used (i.e., so cells with \code{NA} will not have points located within them).
#' @param tol1 Numeric >0, maximum root-mean-square distance allowed between the set of observed pairwise distances between points in \code{x1} and the set of randomized pairwise distances between points simulating \code{x1}. The algorithm will shuffle points until the calculated difference is <= this number. Units are the same as units used by the coordinate reference system of \code{x} (usually meters).
#' @param tol12 As \code{tol1} but for the root-mean-square deviation between points in \code{x1} and \code{x2}.
#' @param distFunct Either a function or \code{NULL}. If \code{NULL} then \code{\link[geosphere]{distGeo}} is used to calculate distances.  Other "dist" functions (e.g., \code{\link[geosphere]{distGeo}}) can be used.  Alternatively, a custom function can be used so long as its first argument is a 2-column numeric matrix with one row for the x- and y-coordinates of a single point and its second argument is a two-column numeric matrix with one or more rows of other points.
#' @param keepData Logical, if \code{TRUE} then the original data in \code{x} (i.e., columns that do not represent coordinates) will be retained in the output but the coordinates will be shuffled. If \code{FALSE} (default) then the returned value will have just shuffled coordinates.
#' @param verbose Logical, if \code{FALSE} (default) show no progress indicator. If \code{TRUE} then display updates and graph.
#' @param ... Arguments to pass to \code{distGeo} or \code{\link[enmSdm]{sampleRast}}. Note that if \code{x} is a matrix or data frame a coordinate reference system may be passed using \code{crs = <proj4 string code>} or \code{crs = <object of class CRS>} (see \pkg{sp} package). Otherwise the coordinates are assumed to be unprojected (WGS84).
#' @return Object of the same class as \code{x1} but with coordinates randomized.
#' @seealso \code{\link[enmSdm]{randPointsRespectingSelf}}, \code{\link[enmSdm]{randPointsRespectingSelfOther2}}, \code{\link[enmSdm]{randPointsBatch}}, \code{\link[enmSdm]{randPointsBatchSampled}}, \code{\link[enmSdm]{randPointsBatchExtract}}, \code{\link[enmSdm]{randPointsBatchNicheOverlap}}
#' @examples
#' library(dismo)
#' library(raster)
#' 
#' data(lemurs, package='enmSdm')
#' longLat <- c('decimalLongitude', 'decimalLatitude')
#' 
#' mad- <- raster::getData('GADM', country='MDG', level=0)
#' elev <- raster::getData('alt', country='MDG', mask=TRUE, res=2.5)
#' 
#' # plot data as-is
#' plot(mad)
#' species <- sort(unique(lemurs$species))
#' 
#' for (i in seq_along(species)) {
#' 	
#' 	thisLemur <- lemurs[lemurs$species == species[i], longLat]
#' 	points(thisLemur, pch=i, col=i)
#' 	
#' }
#' 
#' legend('bottomleft', legend=species, pch=seq_along(species), col=seq_along(species))
#' 
#' # geographically thin presences of each species
#' thinLemurs <- data.frame()
#' 
#' for (i in seq_along(species)) {
#' 	
#' 	thisLemur <- lemurs[lemurs$species == species[i], ]
#' 	thinned <- geoThin(thisLemur, minDist=10000, longLat=longLat)
#' 	thinLemurs <- rbind(thinLemurs, thinned)
#' 	
#' }
#' 
#' # plot geographically thinned data
#' plot(mad)
#' 
#' for (i in seq_along(species)) {
#' 	
#' 	thisLemur <- thinLemurs[thinLemurs$species == species[i], longLat]
#' 	points(thisLemur, pch=i, col=i)
#' 	
#' }
#' 
#' legend('bottomleft', legend=species, pch=seq_along(species), col=seq_along(species))
#' 
#' # randomize one species with respect to itself 
#' x <- thinLemurs[thinLemurs$species == 'Eulemur fulvus', longLat]
#' 
#' set.seed(123)
#' x1rand <- randPointsRespectingSelf(x=x, rast=elev, tol=24000, verbose=TRUE)
#' 
#' # plot observed and randomized occurrences
#' plot(mad)
#' points(x, pch=16)
#' points(x1rand, col='red')
#' 
#' # randomize two species with respect to selves and others
#' species1 <- species[1]
#' species2 <- species[3]
#' 
#' x1 <- thinLemurs[thinLemurs$species == species1, longLat]
#' x2 <- thinLemurs[thinLemurs$species == species2, longLat]
#' 
#' set.seed(123)
#' tol1 <- tol2 <- tol12 <- 16000
#' x12rand <- randPointsRespectingSelfOther2(x1=x1, x2=x2, rast=elev,
#' 	tol1=tol1, tol2=tol2, tol12=tol12, verbose=TRUE)
#' 
#' # plot geographically thinned data
#' plot(mad)
#' points(x1, pch=21, bg='cornflowerblue')
#' points(x2, pch=24, bg='cornflowerblue')
#' points(x12rand$x1rand, pch=1, col='red')
#' points(x12rand$x2rand, pch=2, col='red')
#' 
#' legend('bottomleft', legend=c(species1, species2,
#' 	legend=paste('rand', species1), paste('rand', species2)),
#' 	pch=c(21, 24, 1, 2), col=c('black', 'black', 'red', 'red'),
#' 	pt.bg=c('cornflowerblue', 'cornflowerblue', NA, NA))
#' 
#' ### batch mode
#' \donttest{
#' 
#' # download climate data
#' clim <- raster::getData('worldclim', var='bio', res=2.5)
#' 
#' # lemur data
#' data(lemurs, package='enmSdm')
#' longLat <- c('decimalLongitude', 'decimalLatitude')
#' 
#' # geographically thin presences of each species
#' thinLemurs <- data.frame()
#' 
#' for (i in seq_along(species)) {
#' 	
#' 	thisLemur <- lemurs[lemurs$species == species[i], ]
#' 	thinned <- geoThin(thisLemur, minDist=10000, longLat=longLat)
#' 	thinLemurs <- rbind(thinLemurs, thinned)
#' 	
#' }
#' 
#' # randomize two species with respect to selves and others
#' species1 <- species[1]
#' species2 <- species[3]
#' 
#' x1 <- thinLemurs[thinLemurs$species == species1, longLat]
#' x2 <- thinLemurs[thinLemurs$species == species2, longLat]
#' 
#' # create null distributions
#' set.seed(123)
#' tol1 <- tol2 <- tol12 <- 24000
#' iterations <- 100 # for analysis set this to 100 or more
#' # for testing use a small number!
#' 
#' x12rand <- randPointsBatch('randPointsRespectingSelfOther2', x1=x1, x2=x2,
#' 	rast=clim[[1]], tol1=tol1, tol2=tol2, tol12=tol12, iterations=iterations,
#' 	verbose=TRUE)
#' 
#' # get environment that was sampled to use as background
#' bg <- randPointsBatchSampled(x12rand)
#' bgEnv <- raster::extract(clim, bg)
#' 
#' # create PCA of environmental space
#' vars <- paste0('bio', 1:19)
#' bgPca <- princomp(bgEnv[ , vars], cor=TRUE)
#' 
#' x1env <- raster::extract(clim, x1)
#' x2env <- raster::extract(clim, x2)
#' 
#' nas1 <- omnibus::naRows(x1env)
#' nas2 <- omnibus::naRows(x2env)
#' 
#' if (length(nas1) > 0) x1env <- x1env[-nas1, ]
#' if (length(nas2) > 0) x2env <- x2env[-nas2, ]
#' 
#' # observed niche overlap
#' obsOverlap <- enmSdm::nicheOverlap(
#' 	x1=x1env,
#' 	x2=x2env,
#' 	env=bgPca,
#' 	vars=vars,
#' 	bins=100,
#' 	cor=TRUE
#' )
#' 
#' # extract climate at randomized sites
#' x12rand <- randPointsBatchExtract(x12rand, clim, verbose=TRUE)
#' 
#' # null niche overlap
#' nullOverlap <- randPointsBatchNicheOverlap(
#' 	rands=x12rand,
#' 	env=bgPca,
#' 	vars=vars,
#' 	bins=100,
#' 	cor=TRUE
#' )
#' 
#' hist(nullOverlap$d, 20, main='Niche Overlap',
#' 	xlab='Schoener\'s D', xlim=c(0, 1))
#' abline(v=obsOverlap[['d']], col='blue', lwd=3)
#' legend('topright', legend='Observed', lwd=3, col='blue')
#' 
#' }
#' @export

randPointsRespectingSelfOther1 <- function(
	x1,
	x2,
	rast,
	tol1 = NULL,
	tol12 = NULL,
	distFunct = NULL,
	keepData = FALSE,
	verbose=FALSE,
	...
) {

	ellipses <- list(...)

	if (is.null(distFunct)) distFunct <- geosphere::distGeo

	# save copy of original
	out1 <- x1

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

	# remember CRS
	crs <- if ('crs' %in% ellipses) {
		ellipses$crs
	} else {
		raster::projection(x1)
	}

	### calculate observed pairwise distances
	#########################################

	obsDists1 <- geosphere::distm(x1, fun=distFunct)
	obsDists12 <- geosphere::distm(x1, x2, fun=distFunct)
	
	obsDists1[upper.tri(obsDists1, diag=TRUE)] <- NA

	calcTol <- function(dists, removeTopRow) {
	
		# removeTopRow	Logical, if TRUE remove top row
		# 				(useful to avoid warning if top row
		#				is all NAs)
		
		indices <- if (removeTopRow) {
			2:nrow(dists)
		} else {
			1:nrow(dists)
		}
		
		minDists <- apply(dists[indices, ], 1, min, na.rm=TRUE)
		minPosDists <- minDists[minDists > 0]
		tol <- min(minPosDists)
		tol
	}
	
	if (is.null(tol1)) {
	
		if (verbose) omnibus::say('Calculating tolerance "tol1" automatically using mean of\n minimum pairwise distances for all distances > 0.')
		tol1 <- calcTol(obsDists1, TRUE)
		if (verbose) omnibus::say('Using ', sprintf('%.2f', tol1), ' for tolerance "tol1".', post=2)
	
	}

	if (is.null(tol12)) {
	
		if (verbose) omnibus::say('Calculating tolerance "tol12" automatically using mean of\n minimum pairwise distances for all distances > 0.')
		tol12 <- calcTol(obsDists12, FALSE)
		if (verbose) omnibus::say('Using ', sprintf('%.2f', tol12), ' for tolerance "tol12".', post=2)
	
	}

	### initiate random points
	##########################
	
	x1size <- length(x1)
	x2size <- length(x2)

	x1index <- seq_along(x1)
	
	rastSize <- raster::ncell(rast)
	numRandPoints1 <- min(round(rastSize / 10), max(10000, round(length(x1)^2 / max(1, log(tol1)))))
	numRandPoints2 <- min(round(rastSize / 10), max(10000, round(length(x2)^2 / max(1, log(tol1)))))
	numRandPoints <- 10 * (numRandPoints1 + numRandPoints2)
	
	randPoints <- enmSdm::sampleRast(rast, numRandPoints, prob=FALSE)
	randPoints <- sp::SpatialPoints(randPoints, sp::CRS(crs))
	
	randSites1 <- randPoints[seq_along(x1)]
	randUsed <- x1size

	randDists1 <- geosphere::distm(randSites1, fun=distFunct)
	randDists12 <- geosphere::distm(randSites1, x2, fun=distFunct)

	randDists1[upper.tri(randDists1, diag=TRUE)] <- NA

	delta1 <- statisfactory::rmsd(obsDists1, randDists1, na.rm=TRUE)
	delta12 <- statisfactory::rmsd(obsDists12, randDists12, na.rm=TRUE)
	
	# if (verbose) {

		# histBreaks <- pretty(c(0, max(randDists, obsDists, na.rm=TRUE)), 20)
		# randHist <- hist(randDists, breaks=histBreaks, plot=FALSE)
		# obsHist <- hist(obsDists, breaks=histBreaks, plot=FALSE)
		# plot(obsHist, xlab='Distance', main='Pairwise Distance Distribution')
		# mids <- randHist$mids
		# points(mids, randHist$counts, col='red', pch=5)
		# legend('topright', inset=0.01, legend=c('Observed', 'Current randomized', 'Past randomized'), col=c(NA, 'red', 'gray'), pch=c(NA, 1, 1), fill=c('white', NA, NA), border=c('black', NA, NA))

	# }

	# get format styling for reporting
	if (verbose) {
			
		diff1 <- round(delta1 - tol1)
		diff12 <- round(delta12 - tol12)
	
		ncharDiffs <- nchar(c(as.character(diff1), as.character(diff12)))
		ncharDiffs <- max(ncharDiffs)

		style <- paste0('%', ncharDiffs, '.2f')
				
	}
		
	tries <- accepts <- 0
	
	# iteratively randomized, check to see if this made distribution of randomized distances closer to observed distribution, if so keep
	while ((delta1 > tol1) | (delta12 > tol12)) {

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

		candRandDists <- randDists1
		candRandSites <- randSites1

		replaceIndex <- sample(x1index, 1)
		candDist <- geosphere::distm(candSite, randSites1, fun=distFunct)
		
		candRandDists[replaceIndex, ] <- candRandDists[ , replaceIndex] <- candDist
		candRandDists[upper.tri(candRandDists, diag=TRUE)] <- NA
		candRandSites@coords[replaceIndex, ] <- sp::coordinates(candSite)

		candDists12 <- geosphere::distm(candRandSites, x2, fun=distFunct)
		
		candDelta1 <- statisfactory::rmsd(obsDists1, candRandDists, na.rm=TRUE)
		candDelta12 <- statisfactory::rmsd(obsDists12, candDists12, na.rm=TRUE)

		### accept randomized point
		if ((candDelta1 <= delta1) & (candDelta12 <= delta12)) {

			randSites1@coords[replaceIndex, ] <- sp::coordinates(candSite)
			randDists1 <- randDists1
			randDists12 <- randDists12

			delta1 <- candDelta1
			delta12 <- candDelta12
			
			accepts <- accepts + 1

			# report
			if (verbose) {
			
				diff1 <- delta1 - tol1
				diff12 <- delta12 - tol12
				
				omnibus::say('(actual - desired) tolerances for tol1 and tol12: ', sprintf(style, diff1), ' | ', sprintf(style, diff12))
				
			}
			
		}

	}

	out <- .randPointsReturn(x=x1, randSites=randSites1, crs=crs, keepData=keepData, out=out)
	
	attr(out, 'keepData') <- keepData
	attr(out, 'tries') <- tries
	attr(out, 'accepts') <- accepts
	
	tols <- c(tol1, tol12)
	names(tols) <- c('tol1', 'tol12')
	attr(out, 'tol') <- tols
	
	finalTol <- c(delta1, delta12)
	names(finalTol) <- c('tol1', 'tol12')
	attr(out, 'finalTol') <- finalTol
	
	out

}
