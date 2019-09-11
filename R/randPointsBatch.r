#' Iterate "randPointsRespecting..." functions.
#'
#' This function is a wrapper function for any of the \code{randPointsRespecting~} functions. It is useful for calling the function multiple times using the same arguments. The output is a list with one element per call.
#' @param randFunct Character, any of: \code{'randPointsRespectingSelf'}, \code{'randPointsRespectingSelfOther1'}, or \code{'randPointsRespectingSelfOther2'}.
#' @param ... Arguments to pass to \code{randPointsRespectingSelf}, \code{randPointsRespectingSelfOther1}, or \code{randPointsRespectingSelfOther2}.
#' @param iterations Positive integer, number of times to call the function.
#' @param rast Raster, RasterStack, or RasterBrick used to locate presences randomly. If this is a RasterStack or a RasterBrick then the first layer will be used (i.e., so cells with \code{NA} will not have points located within them).
#' @param distFunct Either a function or \code{NULL}. If \code{NULL} then \code{\link[geosphere]{distGeo}} is used to calculate distances.  Other "dist" functions (e.g., \code{\link[geosphere]{distGeo}}) can be used.  Alternatively, a custom function can be used so long as its first argument is a 2-column numeric matrix with one row for the x- and y-coordinates of a single point and its second argument is a two-column numeric matrix with one or more rows of other points.
#' @param keepData Logical, if \code{TRUE} then the original data in \code{x} (i.e., columns that do not represent coordinates) will be retained in the output but the coordinates will be shuffled. If \code{FALSE} (default) then the returned value will have just shuffled coordinates.
#' @param verbose Logical, if \code{TRUE} (default) show progress along iterations.
#' @param verboseEach Logical, if \code{FALSE} (default), then suppress progress indicators for each call of the \code{randPointsRespecting~} function.
#' @details Note that if you use the \code{randPointsRespectingSelfOther2} function and intend to summarize across iterations using subsequent functions (e.g., \code{\link[enmSdm]{randPointsSampled}}), it is highly advisable to ensure that the \code{x1} and \code{x2} arguments to that function have the same class (matrix, data frame, SpatialPoints, or SpatialPointsDataFrame). You should also either set the argument \code{keepData} to \code{FALSE} or ensure that the column names of \code{x1} and \code{x2} are exactly the same (which then allows you to use \code{keepData = TRUE}).
#' @return A list with \code{iterations} elements.
#' @seealso \code{\link[enmSdm]{randPointsRespectingSelf}}, \code{\link[enmSdm]{randPointsRespectingSelfOther1}}, \code{\link[enmSdm]{randPointsRespectingSelfOther2}}, \code{\link[enmSdm]{randPointsBatchSampled}}, \code{\link[enmSdm]{randPointsBatchExtract}}, \code{\link[enmSdm]{randPointsBatchNicheOverlap}}
#' @examples
#' library(dismo)
#' library(raster)
#'
#' data(lemurs, package='enmSdm')
#' longLat <- c('decimalLongitude', 'decimalLatitude')
#'
#' mad <- raster::getData('GADM', country='MDG', level=0)
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
#' }
#' @export

randPointsBatch <- function(
	randFunct,
	iterations = 100,
	...,
	rast,
	distFunct = NULL,
	keepData = FALSE,
	verbose=TRUE,
	verboseEach=FALSE
) {

	if (verboseEach) verbose <- TRUE

	if (verbose) omnibus::say('Executing ', iterations, ' iterations of ', randFunct, '...')

	out <- list()

	funct <- match.fun(randFunct)

	# for each iteration
	for (iter in 1:iterations) {

		if (verbose) {
			level <- if (verboseEach) { 1} else { NULL }
			omnibus::say('Iteration ', iter, ' of ', iterations, level=level)
		}

		out[[iter]] <- funct(..., rast = rast, distFunct = distFunct, verbose = verboseEach)

	} # next iteration

	attr(out, 'randFunct') <- funct
	attr(out, 'randFunctName') <- randFunct
	attr(out, 'iterations') <- iterations

	out

}
