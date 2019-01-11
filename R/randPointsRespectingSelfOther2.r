#' Randomizes the location of two sets of geographic points while respecting spatial autocorrelation
#'
#' This function randomizes the location of two sets of geographic points with respect to one another retaining (more or less) the same distribution of pairwise distances between points with and between sets (plus or minus a user-defined tolerance).
#' @param x1 Matrix, data frame, SpatialPoints, or SpatialPointsDataFrame object. If this is a matrix or data frame, the first two columns must represent longitude and latitude (in that order). If \code{x} is a matrix or data frame, the coordinates are assumed to be unprojected (WGS84) (a coordinate reference system proj4 string or \code{CRS} object can be passed into the function using \code{...}). If \code{x} is a SpatialPoints or SpatialPointsDataFrame and not in WGS84 or NAD83, then coordinates are projected to WGS84 (with a warning).
#' @param x2 As \code{x1}.
#' @param rast Raster, RasterStack, or RasterBrick used to locate presences randomly. If this is a RasterStack or a RasterBrick then the first layer will be used (i.e., so cells with \code{NA} will not have points located within them).
#' @param tol1 Numeric >0, maximum root-mean-square distance allowed between the set of observed pairwise distances between points in \code{x1} and the set of randomized pairwise distances between points simulating \code{x1}. The algorithm will shuffle points until the calculated difference is <= this number. Units are the same as units used by the coordinate reference system of \code{x} (usually meters).
#' @param tol2 As \code{tol1} but for \code{x2}.
#' @param tol12 As \code{tol1} but for the root-mean-square deviation between points in \code{x1} and \code{x2}.
#' @param distFunct Either a function or \code{NULL}. If \code{NULL} then \code{\link[geosphere]{distGeo}} is used to calculate distances.  Other "dist" functions (e.g., \code{\link[geosphere]{distGeo}}) can be used.  Alternatively, a custom function can be used so long as its first argument is a 2-column numeric matrix with one row for the x- and y-coordinates of a single point and its second argument is a two-column numeric matrix with one or more rows of other points.
#' @param keepData Logical, if \code{TRUE} then the original data in \code{x1} and \code{x1} (i.e., columns that do not represent coordinates) will be retained in the output but the coordinates will be shuffled. If \code{FALSE} (default) then the returned values will have just shuffled coordinates.
#' @param verbose Logical, if \code{FALSE} (default) show no progress indicator. If \code{TRUE} then display updates and graph.
#' @param ... Arguments to pass to \code{distGeo} or \code{\link[enmSdm]{sampleRast}}. Note that if \code{x} is a matrix or data frame a coordinate reference system may be passed using \code{crs = <proj4 string code>} or \code{crs = <object of class CRS>} (see \pkg{sp} package). Otherwise the coordinates are assumed to be unprojected (WGS84).
#' @return A list with two elements, each representing object of the same classes as \code{x1} and \code{x2} but with coordinates randomized.
#' @seealso \code{\link[enmSdm]{randPointsRespectingSelf}}, \code{\link[enmSdm]{randPointsRespectingSelfOther1}}, \code{\link[enmSdm]{randPointsBatch}}, \code{\link[enmSdm]{randPointsBatchSampled}}, \code{\link[enmSdm]{randPointsBatchExtract}}, \code{\link[enmSdm]{randPointsBatchNicheOverlap}}
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

randPointsRespectingSelfOther2 <- function(
	x1,
	x2,
	rast,
	tol1 = NULL,
	tol2 = NULL,
	tol12 = NULL,
	distFunct=NULL,
	restrict = TRUE,
	keepData = FALSE,
	verbose=FALSE,
	...
) {

	ellipses <- list(...)

	if (is.null(distFunct)) distFunct <- geosphere::distGeo

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

	# remember CRS
	crs <- if ('crs' %in% ellipses) {
		ellipses$crs
	} else {
		raster::projection(x1)
	}

	### calculate observed pairwise distances
	#########################################

	obsDists1 <- geosphere::distm(x1, fun=distFunct)
	obsDists2 <- geosphere::distm(x2, fun=distFunct)
	obsDists12 <- geosphere::distm(x1, x2, fun=distFunct)
	
	obsDists1[upper.tri(obsDists1, diag=TRUE)] <- NA
	obsDists2[upper.tri(obsDists2, diag=TRUE)] <- NA

	if (is.null(tol1)) {
	
		tol1 <- .calcTol(obsDists1, TRUE)
		if (verbose) omnibus::say('Automatically calculated value for tol1:  ', sprintf('%.2f', tol1), '.')
	
	}

	if (is.null(tol2)) {
	
		tol2 <- .calcTol(obsDists2, TRUE)
		if (verbose) omnibus::say('Automatically calculated value for tol2:  ', sprintf('%.2f', tol2), '.')
	
	}

	if (is.null(tol12)) {
	
		tol12 <- .calcTol(obsDists12, FALSE)
		if (verbose) omnibus::say('Automatically calculated value for tol12: ', sprintf('%.2f', tol12), '.')
	
	}

	### initiate random points
	##########################
	
	x1size <- length(x1)
	x2size <- length(x2)
	x12size <- x1size + x2size

	x1index <- seq_along(x1)
	x2index <- seq_along(x2)
	
	rastSize <- raster::ncell(rast)
	numRandPoints <- 5 * max(10000, round((length(x1) * length(x2)) / max(1, log(min(tol1, tol2, tol12)))))
	
	x1x2 <- rbind(x1, x2)
	tol <- min(tol1, tol2, tol12)
	
	if (verbose) omnibus::say('Placing points in a spatially structured manner using minimum of tol1, tol2, and tol12...')
	
	randPoints12 <- randPointsRespectingSelf(
		x = x1x2,
		rast = rast,
		tol = tol,
		distFunct = distFunct,
		restrict = restrict,
		keepData = keepData,
		verbose = verbose,
		...
	)
	
	### calculate summary statistics for these point sets
	#####################################################
	
	randSites1 <- randPoints12[x1index]
	randSites2 <- randPoints12[x1size + x2index]
	
	randDists1 <- geosphere::distm(randSites1, fun=distFunct)
	randDists2 <- geosphere::distm(randSites2, fun=distFunct)
	randDists12 <- geosphere::distm(randSites1, randSites2, fun=distFunct)

	randDists1[upper.tri(randDists1, diag=TRUE)] <- NA
	randDists2[upper.tri(randDists2, diag=TRUE)] <- NA

	delta1 <- statisfactory::rmsd(obsDists1, randDists1, na.rm=TRUE)
	delta2 <- statisfactory::rmsd(obsDists2, randDists2, na.rm=TRUE)
	delta12 <- statisfactory::rmsd(obsDists12, randDists12, na.rm=TRUE)

	# tries <- accepts <- 0
	
	# ### iteratively randomize
	# #########################

	# # swap randomly selected sites between species to see if this increases
	# # match to observed summary statistics
	
	# if (verbose) omnibus::say('Swapping assignments between species...', pre=2)
	
	# while ((delta1 > tol1) | (delta2 > tol2) | (delta12 > tol12)) {

		# tries <- tries + 1
		
		# # randomly choose one point from each set and swap
		# # if this reduced all deltas, then accept
		
		# allCandSites1 <- randSites1
		# allCandSites2 <- randSites2
		
		# which1 <- sample(x1index, 1)
		# which2 <- sample(x2index, 1)
		
		# candSite1 <- randSites1[which1]
		# candSite2 <- randSites2[which2]
		
		# allCandSites1@coords[which1, ] <- sp::coordinates(candSite2)
		# allCandSites2@coords[which2, ] <- sp::coordinates(candSite1)
		
		# candDists1 <- geosphere::distm(allCandSites1, x1, fun=distFunct)
		# candDists2 <- geosphere::distm(allCandSites2, x2, fun=distFunct)
		# candDists12 <- geosphere::distm(allCandSites1, allCandSites2, fun=distFunct)
		
		# candDelta1 <- statisfactory::rmsd(candDists1, obsDists1)
		# candDelta2 <- statisfactory::rmsd(candDists2, obsDists2)
		# candDelta12 <- statisfactory::rmsd(candDists12, obsDists12)
		
		# ### accept randomized point
		# if ((candDelta1 <= delta1) & (candDelta2 <= delta2) & (candDelta12 <= delta12)) {

			# randSites1 <- candSites1
			# randSites2 <- candSites2
		
			# delta1 <- candDelta1
			# delta2 <- candDelta2
			# delta12 <- candDelta12
			
			# # report
			# if (verbose) {
			
				# diff1 <- delta1 - tol1
				# diff2 <- delta2 - tol2
				# diff12 <- delta12 - tol12
				
				# omnibus::say('(actual - desired) tolerances for tol1, tol2, and tol12: ', sprintf('%.2f', diff1), ' | ', sprintf('%.2f', diff2), ' | ', sprintf('%.2f', diff12), ' with ', accepts, ' accepts of ', tries, ' tries')
				
			# }
			
		# }

	# }

	out1 <- .randPointsReturn(x=x1, randSites=randSites1, crs=crs, keepData=keepData, out=out1) 
	out2 <- .randPointsReturn(x=x2, randSites=randSites2, crs=crs, keepData=keepData, out=out2) 
	
	out <- list(x1rand=out1, x2rand=out2)

	attr(out, 'keepData') <- keepData
	attr(out, 'distFunct') <- distFunct
	attr(out, 'restrict') <- restrict
	
	tols <- c(tol1, tol2, tol12)
	names(tols) <- c('tol1', 'tol2', 'tol12')
	attr(out, 'tol') <- tols
	
	finalTol <- c(delta1, delta2, delta12)
	names(finalTol) <- c('tol1', 'tol2', 'tol12')
	attr(out, 'finalTol') <- finalTol
	
	out

}
