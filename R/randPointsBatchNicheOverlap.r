#' Niche overlap for a set of iterated "randPointsRespecting~" functions
#'
#' This function is called using a list object typically generated using the \code{\link[enmSdm]{randPointsBatchExtract}} function (plus sometimes followed by the \code{\link[enmSdm]{randPointsBatchExtract}} and \code{\link[enmSdm]{randPointsBatchSampled}} functions). It calculates niche overlap for each set of points. Essentially it is a wrapper for \code{\link[enmSdm]{nicheOverlap}}.
#' @param rands A list object typically generated using the \code{\link[enmSdm]{randPointsBatchExtract}} function.
#' @param env Either a data frame, matrix, or any object that can be coerced to a data frame containing environmental data at available background sites, \emph{or} an object of class \code{princomp} representing a principal components analysis generated using the \code{\link[stats]{princomp}} function with argument \code{scores = TRUE}.
#' @param vars Either a character list naming columns in \code{x1}, \code{x2}, and \code{x3} to be used as environmental data, \emph{or} positive integers indexing the columns to be used as environmental data.
#' @param x Either \code{NULL} (default) or a data frame, matrix, SpatialPointsDataFrame, or other object that can be coerced to a data frame. If supplied then the objects (usually species) represented in \code{rands} are compared to this set of environmental values. This argument \emph{must} be supplied if \code{rands} was generated using \code{\link[enmSdm]{randPointsRespectingSelf}} or \code{\link[enmSdm]{randPointsRespectingSelfOther1}}.
#' @param bins Number of bins into which to divide the environmental space (default is 100 on each side).
#' @param cor Logical, if \code{TRUE} (default), then the PCA used to construct the environmental space will use the correlation matrix (this is highly recommended if the variables are on different scales). This is ignored if \code{env} is an object of class \code{princomp}.
#' @return Data frame.
#' @seealso \code{\link[enmSdm]{randPointsRespectingSelf}}, \code{\link[enmSdm]{randPointsRespectingSelfOther1}}, \code{\link[enmSdm]{randPointsRespectingSelfOther2}}, \code{\link[enmSdm]{randPointsBatch}}, \code{\link[enmSdm]{randPointsBatchSampled}}, \code{\link[enmSdm]{randPointsBatchExtract}}
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
#' \dontrun{
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

randPointsBatchNicheOverlap <- function(
	rands,
	env,
	vars,
	x = NULL,
	bins = 100,
	cor = TRUE
) {

	if (attr(rands, 'randFunctName') %in% c('randPointsRespectingSelf', 'randPointsRespectingSelfOther1')) {
		if (is.null(x)) stop('Argument "x" must be specified if argument "rands" was generated using either "randPointsRespectingSelf" or "randPointsRespectingSelfOther1".')
	}

	out <- data.frame()

	for (i in seq_along(rands)) {

		if (attr(rands, 'randFunctName') %in% c('randPointsRespectingSelf', 'randPointsRespectingSelfOther1')) {
			x1 <- rands[[i]]
			x2 <- x
		} else {
			x1 <- rands[[i]]$x1rand
			x2 <- rands[[i]]$x2rand
		}

		thisOverlap <- nicheOverlap(x1=x1, x2=x2, env=env, vars=vars, bins=bins, cor=cor)
		thisOverlap <- as.data.frame(thisOverlap)
		thisOverlap <- t(thisOverlap)
		out <- rbind(out, thisOverlap)

	}

	out <- cbind(
		data.frame(iter = 1:nrow(out)),
		out
	)
	
	rownames(out) <- 1:nrow(out)

	out

}
