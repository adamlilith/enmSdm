#' Local characteristic distance of spatial autocorrelation for variables
#'
#' This function calculates the range of spatial autocorrelation for a set of predictors at a set of sites. Input can be a raster or raster stack, in which case the output is a raster stack, one layer per layer of input, with cell values equal to the characteristic distance of spatial autocorrelation for each cell.  Alternatively, input can be a raster and a matrix, data frame, SpatialPoints, or SpatialPoointsDataFrame object, in which case the output will be a matrix, data frame, or SpatialPointsDataFrame with the characteristic distance of spatial autocorrelation for each layer in the raster set at each point. Finally, input can simply be a matrix, data frame, of SpatialPointsDataFrame in which case the scale of autocorrelation is calculated using data from the sites, with output format matching the input format. See \emph{Details} for information on how the characteristic scale of spatial autocorrelation is estimated. This function is related to \code{\link[enmSdm]{spatialCorrForPoints}} which calculates spatial autocorrelation for distances between points, whereas this function calculates spatial autocorrelation for measurements taken at points (or raster cell centers).
#' @param x Either a raster, raster stack/brick, matrix with column names, data frame, or SpatialPointsDataFrame. If you use a matrix or data frame then the first two columns will be assumed to represent longitude and latitude, in that order, and their coordinate reference system will be assumed to be WGS84 (unprojected).
#' @param focal This has various uses depending on the type of object specified for \code{x}:
#' \itemize{
#'		\item If \code{x} is a raster or raster stack/brick and \code{focal} is \code{NULL} (default), then spatial autocorrelation will be calculated across all non-\code{NA} cells in  \code{x}.
#'		\item If \code{x} is a SpatialPointsDataFrame and \code{focal} is \code{NULL} (default), then spatial autocorrelation will be calculated across all non-\code{NA} points represented by \code{x}.
#' 		\item If \code{x} is a raster or raster stack/brick and \code{x} is a matrix with column names, data frame, SpatialPoints, or SpatialPointsDataFrame object, then autocorrelation will be calculated focal all sites in \code{focal}. If you use a matrix or data frame then the first two columns will be assumed to represent longitude and latitude, in that order, and their coordinate reference system will be assumed to be the same as the raster or raster stack/brick.
#' }
#' @param vars Character vector or \code{NULL} (default). If \code{x} is a SpatialPointsDataFrame, by default the variables on which to calculate spatial autocorrelation will be taken from the column names of \code{x}. However, you can use \code{vars} to specify the column names in \code{x} to use. This is ignored if \code{x} is not a SpatialPointsDataFrame.
#' @param breaks Either: A single integer, three numeric values, or a matrix or data frame with at least two columns:
#' \itemize{
#' 		\item A positive integer: Number of overlapping bins to use. The minimum bin distance will be >0 and maximum <= the maximum inter-point distance.
#' 		\item Three numeric values: The first two values are smallest and largest distances (in units used in coordinate reference system of \code{rast}, typically meters) across which to tabulate distance frequencies. The third value is the number of bins.
#' 		\item Matrix or data frame with at least two columns. Each row corresponds to a different bin. The first column represents the minimum distance in each bin (in units used in coordinate reference system of \code{rast}, typically meters) and the second column the maximum distance. Subsequent columns are ignored. Note that by using this option arbitrary bins can be used--they need not overlap or even be continuous in coverage.
#' }
#' @param limitDist Logical, if \code{TRUE}, consider only inter-point distances less than the maximum distance indicated by \code{breaks} (if \code{breaks} is a matrix, data frame, or a 3-element numeric vector). If \code{FALSE} (default), consider all inter-point distances even if they fall beyond the maximum distance defined by  \code{breaks}. Note that if \code{limitDist} is \code{TRUE}, then there is the possibility that a particular cell/point will have no or very few "neighbors". If there are none, the value assigned to the neighbor-less point will be \code{NA}, and if there are few neighbors, then statistical power will be diminished and it will be likely that the maximum distance will be assigned to the cell/point.
#' @param returnMax Logical, if \code{TRUE} (default) and the local characteristic distance of spatial autocorrelation is greater than the maximum distance under consideration, then return the maximum distance (i.e., the upper limit of the last distance bin generated from \code{breaks}). If \code{FALSE} then return \code{NA} in this case.
#' @param iters Positive integer, number of times to permute the randomization test. Default is \code{NULL}.
#' @param perc Numeric value in the range [0, 100], indicates the upper quantile of randomized distance frequencies above which observed distance frequencies are considered significant. Default is 95.
#' @param verbose Logical, if \code{TRUE} (default) show progress.
#' @param ... Other arguments (not used).
#' @return A raster, raster stack/brick, or SpatialPointsDataFrame, depending on the input. If the output is a raster, raster stack, or raster brick, then values represent the minimum distance at which spatial autocorrelation is not significantly different from random expectation (i.e., the "left" side of the distance bin at which this occurs). If the output is a SpatialPointsDataFrame, then columns will be named "sacDistMin_XYZ", "sacDistMid_XYZ" and "sacDistMax_XYZ", where "XYZ" refers to the variable name and "Min", "Mid", and "Max" refer to the minimum, middle, and maximum distance of spatial autocorrelation for each site in the output (i.e., the left, middle, and right side of teh distance bin at which spatial autocorrelation for the site is not different from random expectation).
#' @details The characteristic scale of spatial autocorrelation for a variable for a specific "focal" site relative to a set of "reference" sites is estimated through a multi-step process. The nature of the focal and reference sites depends on the values of \code{x} and \code{focal}. If \code{x} is supplied but \code{focal} is not, then all sites (or cells) in \code{x} will be assumed to be the reference and focal sites. If \code{x} and \code{focal} are supplied, then sites in \code{focal} are the focal sites and sites (or cells) in \code{x} are the reference sites.
#' \enumerate{
#'  \item For each distance bin, calculate the observed upper quantile of the distribution of absolute difference between the value of the variable at the the focal site and all other "reference" sites supplied in argument \code{x}. Here, the upper quantile is given by 100 - 0.5 * (100 - \code{perc}). So if \code{perc = 95} (default), this is the 97.5th quantile of the observed absolute difference for values associated with all reference points in each distance bin.
#'	\item Apply a permutation test by scrambling the absolute differences associated with each pairwise distance between the focal site at reference sites. For each distance bin generate a null expectation from the randomized absolute differences by calculating the value of the 100 - \code{perc}th quantile (so if \code{perc = 95} this is the 5th quantile of the distribution). Repeat \code{iters} times and for each bin calculate the mean of the 100 - \code{perc}th quantile.
#	\item Starting at the first distance bin (smallest distance from the focal site), find the bin at which the observed absolute upper quantile (values from step 1) first falls above the mean of the 100 - \code{perc}th quantile of the randomized values (values from step 2). The distances (lower/middle/upper) associated with the bin at which this occurs represent the characteristic scale of spatial autocorrelation for the given variable at the focal site.
#' }
#' Note that this measure of spatial autocorrelation assumes anisotropy, meaning that from a given focal site the characteristic distance of spatial autocorrelation is the same in all directions.
#' @seealso \code{\link[enmSdm]{spatialCorrForPoints}}
#' @examples
#' \dontrun{
#' # get rasters for mean annual temperature and total precipitation
#' worldClim <- raster::getData('worldclim', var='bio', res=10)
#' worldClim <- raster::subset(worldClim, c(1, 12))
#' 
#' # crop to Madagascar
#' data(mad0)
#' madClim <- raster::crop(worldClim, mad0)
#' 
#' # remove non-Malagasy islands
#' madRast <- raster::rasterize(mad0, madClim)
#' madClim <- madClim * madRast
#' names(madClim) <- c('bio1', 'bio12')
#' 
#' ### spatial autocorrelation for raster (can take a long time!)
#' sacRast <- localSpatialCorrForValues(x=madClim, focal=NULL)
#' sacRast <- sacRast / 1000 # convert to km
#' 
#' # plot
#' par(mfrow=c(2, 2))
#' raster::plot(madClim[['bio1']], main='BIO 01\n(10 * deg C)')
#' raster::plot(madClim[['bio12']], main='BIO 12\n(mm)')
#' raster::plot(sacRast[['bio1']], main='BIO 01\nDistance (km)')
#' raster::plot(sacRast[['bio12']], main='BIO 12\nDistance (km)')
#' 
#' ### spatial autocorrelation for spatial points
#' # faster than others, but still takes a while
#' set.seed(123)
#' x <- dismo::randomPoints(madClim, 200)
#' env <- raster::extract(madClim, x)
#' x <- cbind(x, env)
#' breaks <- c(0, 100000, 10)
#' sacPoints <- localSpatialCorrForValues(x=x, breaks=breaks)
#' 
#' # plot: code point color by characteristic distance of spatial autocorrelation
#' maxSacBio1 <- max(sacPoints$sacDistMid_bio1, na.rm=TRUE)
#' maxSacBio12 <- max(sacPoints$sacDistMid_bio12, na.rm=TRUE)
#' sacBio1 <- 4 * sacPoints$sacDistMid_bio1 / maxSacBio1
#' sacBio12 <- 4 * sacPoints$sacDistMid_bio12 / maxSacBio12
#' 
#' par(mfrow=c(1, 2))
#' leg <- '\n(small: short dist, large: long dist)'
#' sp::plot(mad0, main=paste0('BIO 01', leg))
#' points(sacPoints, pch=1, cex=sacBio1, col='red')
#' sp::plot(mad0, main=paste0('BIO 12', leg))
#' points(sacPoints, pch=1, cex=sacBio12, col='blue')
#'
#' par(mfrow=c(1, 1))
#' maxDist <- max(c(sacPoints$sacDistMid_bio1, sacPoints$sacDistMid_bio12),
#' 	na.rm=TRUE)
#' hist(sacPoints$sacDistMid_bio1, col='red', xlim=c(0, maxDist), breaks=12)
#' hist(sacPoints$sacDistMid_bio12, border='blue', breaks=12, add=TRUE)
#' 
#' ### spatial autocorrelation for spatial points using raster as reference
#' # can take a long time!
#' data(lemur)
#' sacLemur <- localSpatialCorrForValues(x=madClim, focal=lemur)
#' 
#' # plot: code point color by characteristic distance of spatial autocorrelation
#' maxSacBio1 <- max(sacLemur$sacDistMid_bio1, na.rm=TRUE)
#' maxSacBio12 <- max(sacLemur$sacDistMid_bio1, na.rm=TRUE)
#' grayBio1 <- 100 - round(100 * sacLemur$sacDistMid_bio1 / maxSacBio1)
#' grayBio12 <- 100 - round(100 * sacLemur$sacDistMid_bio12 / maxSacBio12)
#' 
#' par(mfrow=c(1, 2))
#' leg <- '\n(dark: short dist, light: long dist)'
#' raster::plot(madClim[['bio1']], cex=1.2, main=paste0('BIO 01', leg))
#' points(lemur, pch=21, bg=paste0('gray', grayBio1))
#' raster::plot(madClim[['bio12']], cex=1.2, main=paste0('BIO 12', leg))
#' points(lemur, pch=21, bg=paste0('gray', grayBio12))
#' 
#' # plot buffers showing characteristic distance of spatial autocorrelation
#' # around each point
#' sacLemurEa <- sp::spTransform(sacLemur, getCRS('mollweide', TRUE))
#' buffsBio1Ea <- rgeos::gBuffer(sacLemurEa, byid=TRUE,
#' 		width=sacLemurEa$sacDistMid_bio1)
#' buffsBio12Ea <- rgeos::gBuffer(sacLemurEa, byid=TRUE,
#' 		width=sacLemurEa$sacDistMid_bio12)
#' buffsBio1 <- sp::spTransform(buffsBio1Ea, getCRS('wgs84', TRUE))
#' buffsBio12 <- sp::spTransform(buffsBio1Ea, getCRS('wgs84', TRUE))
#' par(mfrow=c(1, 2))
#' raster::plot(madClim[['bio1']], main=paste0('BIO 01', leg))
#' sp::plot(buffsBio1, add=TRUE)
#' points(lemur, pch=21, bg=paste0('gray', grayBio1))
#' raster::plot(madClim[['bio12']], main=paste0('BIO 12', leg))
#' sp::plot(buffsBio12, add=TRUE)
#' points(lemur, pch=21, bg=paste0('gray', grayBio12))
#' }
#' @export

localSpatialCorrForValues <- function(
	x,
	focal = NULL,
	vars = NULL,
	breaks = 20,
	limitDist = FALSE,
	returnMax = TRUE,
	iters = 100,
	perc = 95,
	verbose = TRUE,
	...
) {

	focalSupplied <- if (is.null(focal)) { FALSE } else { TRUE }

	# generate common working objects
	# "refs" is set of reference points and associated data
	# "focal" is set of points focal which to calculate SAC plus associated data
	# note these may be the same set of points, depending on the input types
	if (any(c('RasterLayer', 'RasterStack', 'RasterBrick') %in% class(x))) {

		# make a "refs" a SPSDF
		vals <- raster::as.data.frame(x)
		vars <- names(x)
		names(vals) <- vars

		refs <- raster::xyFromCell(x, cell=1:raster::ncell(x), spatial=TRUE)
		refs <- sp::SpatialPointsDataFrame(
			sp::coordinates(refs),
			data=vals,
			proj4=sp::CRS(raster::projection(x))
		)
		
		nonNas <- which(complete.cases(refs@data))
		if (length(nonNas) != raster::ncell(x)) refs <- refs[nonNas, ]

		# convert "focal" to SPSDF
		if (is.null(focal)) {
			
			focal <- refs

		} else if (any(c('matrix', 'data.frame') %in% class(focal))) {
		
			crs <- sp::CRS(raster::projection(x))
		
			focal <- focal[ , 1:2]
			focal <- SpatialPoints(focal, crs)
			vals <- raster::extract(x, focal)
			colnames(vals) <- names(x)
			focal <- sp::SpatialPointsDataFrame(focal, data=vals, crs)
		
		} else if ('SpatialPoints' %in% class(focal)) {
		
			crs <- sp::CRS(raster::projection(focal))

			vals <- raster::extract(x, focal)
			colnames(vals) <- names(x)
			focal <- sp::SpatialPointsDataFrame(focal, data=vals, crs)
			
		} else if ('SpatialPointsDataFrame' %in% class(focal)) {

			vals <- raster::extract(x, focal)
			vals <- as.data.frame(vals)
			colnames(vals) <- names(x)
			focal@data <- vals
			
		}
		
		refsAndFocalSame <- FALSE
		
	} else if (c('SpatialPointsDataFrame') %in% class(x) & is.null(focal)) {

		focal <- refs <- x
		if (is.null(vars)) vars <- names(refs)
		
		refsAndFocalSame <- TRUE
		
	} else if (any(c('matrix', 'data.frame') %in% class(x)) & is.null(focal)) {

		if (is.null(vars)) vars <- colnames(x)[3:ncol(x)]
		data <- as.data.frame(x[ , vars, drop=FALSE])
		x <- sp::SpatialPointsDataFrame(x[ , 1:2], data=data, proj4=enmSdm::getCRS('wgs84', TRUE))
		focal <- refs <- x
		
		refsAndFocalSame <- TRUE
		
	} else {
		stop('Incorrect arguments: "x" must be a raster/raster stack and "focal" NULL; OR "x" must be a raster/raster stack and "focal" a matrix, data frame, SpatialPoints, or SpatialPointsDataFrame; OR "x" must be a matrix, data frame, or SpatialPointsDataFrame, and "focal" NULL.')
	}

	### pre-calculations and -allocations
	
	# number of variables for which to calculate SAC
	numVars <- length(vars)
	
	# pre-calculate quantiles for OBSERVED difference distribution
	twoTailLower <- 0.5 * (100 - perc) / 100
	twoTailMid <- 0.5
	twoTailUpper <- (100 - 0.5 * (100 - perc)) / 100
	
	# pre-calculate quantiles for PERMUTED difference distribution
	oneTailLower <- (100 - perc) / 100
	oneTailMid <- 0.5 * (2 - perc / 100)
	oneTailUpper <- 1
	
	# pre-calculate length of reference vector
	refsSize <- nrow(refs)

	# pre-calculate number of breaks if not supplied
	numBreaks <- if (!any(c('matrix', 'data.frame') %in% class(breaks)) && length(breaks) == 1) {
		breaks
	} else if (!any(c('matrix', 'data.frame') %in% class(breaks)) && length(breaks) == 3) {
		breaks[3]
	} else {
		nrow(breaks)
	}

	# maximum distance to account for
	maxDistToInclude <- if (!any(c('matrix', 'data.frame') %in% class(breaks)) && length(breaks) == 1) {
		Inf
	} else if (!any(c('matrix', 'data.frame') %in% class(breaks)) && length(breaks) == 3) {
		breaks[2]
	} else {
		breaks[numBreaks, 2, drop=TRUE]
	}

	# pre-allocate matrix for permuted lower/middle/upper quantiles of null difference distribution
	randAbsDiffBlank <- matrix(NA, ncol=numBreaks, nrow=iters)
			
	### calculate characteristic scale of spatial autocorrelation for each variable
	###############################################################################
	
	for (countVar in seq_along(vars)) {
	
		thisVar <- vars[countVar]
		
		if (verbose) omnibus::say(thisVar)
			
		# for saving characteristic distance of SAC
		focal@data$DUMMY3 <- focal@data$DUMMY2 <- focal@data$DUMMY1 <- NA
		names(focal@data)[(ncol(focal@data) - 2):ncol(focal@data)] <- paste0(c('sacDistMin_', 'sacDistMid_', 'sacDistMax_'), thisVar)

		# observed values
		focalObs <- focal@data[ , thisVar, drop=TRUE]
		refsObs <- refs@data[ , thisVar, drop=TRUE]

		if (verbose) progress <- txtProgressBar(min=0, max=nrow(focal), style=3, width=min(40, getOption('width')))

		# for each focal point
		for (countFocal in 1:nrow(focal)) {

			# the focal point and its data
			thisFocal <- focal[countFocal, , drop=FALSE]
						
			### distances
			
			# inter-point distances
			dists <- geosphere::distm(thisFocal, refs)
			dists <- c(dists)

			# remove distance to self
			if (refsAndFocalSame) dists[countFocal] <- NA
			
			# limit distances being considered
			if (limitDist) {
				dists <- c(dists[dists <= maxDistToInclude])
				if (length(dists) > 0) dists <- na.omit(dists)
			}

			# if any inter-point distances are in the set being considered
			if (length(dists) > 0) {
				
				# distDistrib <- statisfactory::histOverlap(dists, breaks=breaks, graph=FALSE, indices=TRUE, ...)
				distDistrib <- histOverlap(dists, breaks=breaks, graph=FALSE, indices=TRUE)
				indices <- attr(distDistrib, 'indices')

				# observed differences in values
				obsAbsDiff <- list()
				for (countBin in 1:nrow(distDistrib)) {
					obsAbsDiff[[countBin]] <- if (distDistrib[countBin, 'count'] > 0) {
						abs(refsObs[indices[[countBin]]] - thisFocal@data[ , thisVar])
					} else {
						NA
					}
				}

				# observed lower/middle/upper quantiles of observed difference distribution
				obsAbsDiffLower <- sapply(obsAbsDiff, quantile, p=twoTailLower, na.rm=TRUE)
				# obsAbsDiffMid <- sapply(obsAbsDiff, quantile, twoTailMid, na.rm=TRUE)
				obsAbsDiffMid <- sapply(obsAbsDiff, mean, na.rm=TRUE)
				obsAbsDiffUpper <- sapply(obsAbsDiff, quantile, p=twoTailUpper, na.rm=TRUE)

				# matrices to store lower/middle/upper quantiles of each iterations differences for each distance
				randAbsDiffUpper <- randAbsDiffMid <- randAbsDiffLower <- randAbsDiffBlank
				
				# initialize list to store randomized differences, one per distance bin
				randAbsDiff <- list()

				# for each iteration
				for (iter in 1:iters) {
				
					# scramble reference values
					refsRand <- sample(refsObs, refsSize)
					
					# re-assign scrambled reference values
					for (countBin in 1:nrow(distDistrib)) {
						randAbsDiff[[countBin]] <- if (distDistrib[countBin, 'count'] > 0) {
							abs(refsRand[indices[[countBin]]] - thisFocal@data[ , thisVar])
						} else {
							NA
						}
					}
				
					randAbsDiffLower[iter, ] <- sapply(randAbsDiff, quantile, oneTailLower, na.rm=TRUE)
					# randAbsDiffMid[iter, ] <- sapply(randAbsDiff, quantile, oneTailMid, na.rm=TRUE)
					randAbsDiffMid[iter, ] <- sapply(randAbsDiff, mean, na.rm=TRUE)
					randAbsDiffUpper[iter, ] <- sapply(randAbsDiff, quantile, oneTailUpper, na.rm=TRUE)
				
				} # next iteration
			
				# calculate average lower quantile value for randomized distance distributions
				randAbsDiffLower <- apply(randAbsDiffLower, 2, mean, na.rm=TRUE)
				randAbsDiffMid <- apply(randAbsDiffMid, 2, mean, na.rm=TRUE)
				randAbsDiffUpper <- apply(randAbsDiffUpper, 2, mean, na.rm=TRUE)
		
				### remember the least/middle/most distance where the observed difference in values is no different than the randomized differences
				sacMinDistIndex <- which(randAbsDiffLower <= obsAbsDiffUpper)[1]
				
				# if there is a characteristic distance of SAC
				if (!is.na(sacMinDistIndex)) {
					focal@data[countFocal, paste0('sacDistMin_', thisVar)] <- max(0, distDistrib[sacMinDistIndex, 'lower'])
					focal@data[countFocal, paste0('sacDistMid_', thisVar)] <- max(0, distDistrib[sacMinDistIndex, 'middle'])
					focal@data[countFocal, paste0('sacDistMax_', thisVar)] <- max(0, distDistrib[sacMinDistIndex, 'upper'])
				# if none occurred in the intervals considered, return maximum value (if desired)
				} else if (returnMax) {
					numDistBins <- nrow(distDistrib)
					focal@data[countFocal, paste0('sacDistMin_', thisVar)] <- max(0, distDistrib[numDistBins, 'upper'])
					focal@data[countFocal, paste0('sacDistMid_', thisVar)] <- max(0, distDistrib[numDistBins, 'upper'])
					focal@data[countFocal, paste0('sacDistMax_', thisVar)] <- max(0, distDistrib[numDistBins, 'upper'])
				}

				# ### plot
				# plot(distDistrib[ , 'middle'], obsAbsDiffUpper, ylim=c(0, max(obsAbsDiffUpper,randAbsDiffLower, randAbsDiffUpper, na.rm=TRUE)), lty='dashed')
				
				# polygon(c(distDistrib[ , 'middle'], rev(distDistrib[ , 'middle'])), c(randAbsDiffLower, rev(randAbsDiffUpper)), col='gray')
				# lines(distDistrib[ , 'middle'], randAbsDiffMid, col='black', lwd=2)
				
				# polygon(c(distDistrib[ , 'middle'], rev(distDistrib[ , 'middle'])), c(obsAbsDiffLower, rev(obsAbsDiffUpper)), col=scales::alpha('darkgreen', 0.5))
				# lines(distDistrib[ , 'middle'], obsAbsDiffMid, col='darkgreen', lwd=2)

				# legend('topleft', bty='n', legend=c('random', 'random mean', 'observed', 'observed mean'), fill=c('gray', NA, 'darkgreen', NA, 'darkgreen'), col=c(NA, 'gray', NA, 'darkgreen'), lwd=c(NA, 2, NA, 2))
				
				# abline(v=distDistrib[sacMinDistIndex, 'middle'], col='red')
				
			} # if any inter-point distances in the set being considered
				
			if (verbose) setTxtProgressBar(progress, countFocal)
			
		} # next "focal" site
		
		if (verbose) close(progress)
		
	} # next variable
	
	# put into raster format
	if (any(c('RasterLayer', 'RasterStack', 'RasterBrick') %in% class(x)) & !focalSupplied) {
	
		for (countVar in 1:length(vars)) {
			
			focalVals <- rep(NA, raster::ncell(x))
			if (length(nonNas) > 0) {
				focalVals[nonNas] <- focal@data[ , paste0('sacDistMin_', vars[countVar])]
			} else {
				focalVals <- focal@data[ , paste0('sacDistMin_', vars[countVar])]
			}
		
			thisOut <- x[[1]]
			thisOut <- raster::setValues(thisOut, values=focalVals)
			
			out <- if (exists('out', inherits=FALSE)) {
				raster::stack(out, thisOut)
			} else {
				thisOut
			}
			
		} # next raster layer
			
		if (class(x) == 'RasterBrick') out <- raster::brick(out)
		names(out) <- vars
		
		out

	# else other format
	} else {
	
		focal
		
	}
	
}
