#' Local characteristic distance of spatial autocorrelation for variables
#'
#' This function calculates the range of spatial autocorrelation for a set of predictors at a set of sites. The output is determined by the type of input:
#' \itemize{
#' 		\item If the input is raster or raster stack, the output is a raster stack, one layer per layer of input, with cell values equal to the characteristic distance of spatial autocorrelation for each cell.
#' 		\item If input is a raster \emph{and} a matrix, data frame, \code{SpatialPoints}, or \code{SpatialPointsDataFrame} object, the output will be a matrix, data frame, or \code{SpatialPointsDataFrame} with the characteristic distance of spatial autocorrelation for each layer in the raster set at each point.
#' 		\item If input is just a single matrix, data frame, or \code{SpatialPointsDataFrame}, output will reflect the scale of autocorrelation using data at just the sites represented by the input. The output format be the same as the input format.
#' }
#' See \emph{Details} for information on how the characteristic scale of spatial autocorrelation is estimated. This function is related to \code{\link[enmSdm]{spatialCorrForPoints}} which calculates spatial autocorrelation for distances between points. However, this function calculates spatial autocorrelation for numeric-valued \emph{measurements} taken at points (or raster cell centers).
#' @param x Either a raster, raster stack/brick, matrix with column names, data frame, or SpatialPointsDataFrame. If you use a matrix or data frame then the first two columns will be assumed to represent longitude and latitude, in that order, and their coordinate reference system will be assumed to be WGS84 (unprojected).
#' @param focals This has various uses depending on the type of object specified for \code{x}:
#' \itemize{
#'		\item If \code{x} is a raster or raster stack/brick and \code{focals} is \code{NULL} (default), then spatial autocorrelation will be calculated across all non-\code{NA} cells in  \code{x}.
#'		\item If \code{x} is a SpatialPointsDataFrame and \code{focals} is \code{NULL} (default), then spatial autocorrelation will be calculated across all non-\code{NA} points represented by \code{x}.
#' 		\item If \code{x} is a raster or raster stack/brick and \code{x} is a matrix with column names, data frame, SpatialPoints, or SpatialPointsDataFrame object, then autocorrelation will be calculated focals all sites in \code{focals}. If you use a matrix or data frame then the first two columns will be assumed to represent longitude and latitude, in that order, and their coordinate reference system will be assumed to be the same as the raster or raster stack/brick.
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
#' @details The characteristic scale of spatial autocorrelation for a variable for a specific "focals" site relative to a set of "reference" sites is estimated through a multi-step process. The nature of the focals and reference sites depends on the values of \code{x} and \code{focals}. If \code{x} is supplied but \code{focals} is not, then all sites (or cells) in \code{x} will be assumed to be the reference and focals sites. If \code{x} and \code{focals} are supplied, then sites in \code{focals} are the focals sites and sites (or cells) in \code{x} are the reference sites.
#' \enumerate{
#'  \item For each distance bin, calculate the observed upper quantile of the distribution of absolute difference between the value of the variable at the the focals site and all other "reference" sites supplied in argument \code{x}. Here, the upper quantile is given by 100 - 0.5 * (100 - \code{perc}). So if \code{perc = 95} (default), this is the 97.5th quantile of the observed absolute difference for values associated with all reference points in each distance bin.
#'	\item Apply a permutation test by scrambling the absolute differences associated with each pairwise distance between the focals site at reference sites. For each distance bin generate a null expectation from the randomized absolute differences by calculating the value of the 100 - \code{perc}th quantile (so if \code{perc = 95} this is the 5th quantile of the distribution). Repeat \code{iters} times and for each bin calculate the mean of the 100 - \code{perc}th quantile.
#	\item Starting at the first distance bin (smallest distance from the focals site), find the bin at which the observed absolute upper quantile (values from step 1) first falls above the mean of the 100 - \code{perc}th quantile of the randomized values (values from step 2). The distances (lower/middle/upper) associated with the bin at which this occurs represent the characteristic scale of spatial autocorrelation for the given variable at the focals site.
#' }
#' Note that this measure of spatial autocorrelation assumes anisotropy, meaning that from a given focals site the characteristic distance of spatial autocorrelation is the same in all directions.
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
#' madClim[['bio1']] <- madClim[['bio1']] / 10 # to deg C
#' 
#' ### spatial autocorrelation for raster (can take a long time!)
#' sacRast <- localSpatialCorrForValues(x=madClim, focals=NULL)
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
#' # faster than rasters, but still takes a while
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
#' # plot: notice points along transitional zones
#' # have smaller distances, as expected
#' par(mfrow=c(1, 2))
#' leg <- '\n(small: short dist, large: long dist)'
#' sp::plot(mad0, main=paste0('BIO 01', leg))
#' raster::plot(madClim[['bio1']], add=TRUE)
#' points(sacPoints, pch=1, cex=sacBio1, col='red')
#' sp::plot(mad0, main=paste0('BIO 12', leg))
#' raster::plot(madClim[['bio12']], add=TRUE)
#' points(sacPoints, pch=1, cex=sacBio12, col='blue')
#' #'
#' par(mfrow=c(1, 1))
#' maxDist <- max(c(sacPoints$sacDistMid_bio1, sacPoints$sacDistMid_bio12),
	#' na.rm=TRUE)
#' yMax <- max(
	#' hist(sacPoints$sacDistMid_bio1, breaks=10, plot=FALSE)$counts,
#' 	hist(sacPoints$sacDistMid_bio12, breaks=10, plot=FALSE)$counts
#' )
#' 
#' hist(sacPoints$sacDistMid_bio1, col='red', xlim=c(0, maxDist),
	#' ylim=c(0, yMax), breaks=10, xlab='Distance (m)')
#' hist(sacPoints$sacDistMid_bio12, border='blue', breaks=10, density=10,
	#' add=TRUE)
#' legend('topleft', inset=0.01, legend=c('BIO1', 'BIO12'), fill=c('red', NA),
	#' density=c(NA, 15), border=c('black', 'blue'))
#' 	
#' 
#' ### spatial autocorrelation for spatial points using raster as reference
#' data(lemurs)
#' lemur <- lemurs[lemurs$species == 'Eulemur rufifrons',
	#' c('longitude', 'latitude')]
#' 	
#' # remove record in water
#' lemurBio1 <- raster::extract(madClim[['bio1']], lemur)
#' nas <- which(is.na(lemurBio1))
#' lemur <- lemur[-nas, ]
#' 
#' sacLemur <- localSpatialCorrForValues(x=madClim, focals=lemur, breaks=10)
#' 
#' # plot: code point color by characteristic distance of spatial autocorrelation
#' bio1SacScaled <- 
	#' round(100 * omnibus::stretchMinMax(sacLemur$sacDistMid_bio1))
#' bio12SacScaled <-
	#' round(100 * omnibus::stretchMinMax(sacLemur$sacDistMid_bio12))
#' 
#' grayBio1 <- paste0('gray', bio1SacScaled)
#' grayBio12 <- paste0('gray', bio12SacScaled)
#' 
#' par(mfrow=c(1, 2))
#' leg <- '\n(dark: short dist, light: long dist)'
#' raster::plot(madClim[['bio1']], cex=1.2, main=paste0('BIO 01', leg))
#' points(lemur, pch=21, bg=grayBio1)
#' raster::plot(madClim[['bio12']], cex=1.2, main=paste0('BIO 12', leg))
#' points(lemur, pch=21, bg=grayBio12)
#' 
#' # plot buffers showing characteristic distance of spatial autocorrelation
#' # around each point
#' sacLemurEa <- sp::spTransform(sacLemur, getCRS('mollweide', TRUE))
#' buffsBio1Ea <- rgeos::gBuffer(sacLemurEa, byid=TRUE,
		#' width=sacLemurEa$sacDistMid_bio1)
#' buffsBio12Ea <- rgeos::gBuffer(sacLemurEa, byid=TRUE,
		#' width=sacLemurEa$sacDistMid_bio12)
#' buffsBio1 <- sp::spTransform(buffsBio1Ea, getCRS('wgs84', TRUE))
#' buffsBio12 <- sp::spTransform(buffsBio1Ea, getCRS('wgs84', TRUE))
#' par(mfrow=c(1, 2))
#' raster::plot(madClim[['bio1']], main=paste0('BIO 01', leg))
#' sp::plot(buffsBio1, add=TRUE)
#' points(lemur, pch=21, bg=grayBio1)
#' raster::plot(madClim[['bio12']], main=paste0('BIO 12', leg))
#' sp::plot(buffsBio12, add=TRUE)
#' points(lemur, pch=21, bg=grayBio12)
#' 
#' }

localSpatialCorrForValues <- function(
	x,
	focals = NULL,
	vars = NULL,
	breaks = 20,
	limitDist = FALSE,
	returnMax = TRUE,
	iters = 100,
	perc = 95,
	verbose = TRUE,
	...
) {

	focalSupplied <- if (is.null(focals)) { FALSE } else { TRUE }

	# generate common working objects
	# "refs" is set of reference points and associated data
	# "focals" is set of points at which to calculate SAC plus associated data
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
		
		nonNas <- which(stats::complete.cases(refs@data))
		if (length(nonNas) != raster::ncell(x)) refs <- refs[nonNas, ]

		# convert "focals" to SPSDF
		if (is.null(focals)) {
			
			focals <- refs

		} else if (any(c('matrix', 'data.frame') %in% class(focals))) {
		
			crs <- sp::CRS(raster::projection(x))
		
			focals <- focals[ , 1:2]
			focals <- sp::SpatialPoints(focals, crs)
			vals <- raster::extract(x, focals)
			vals <- as.data.frame(vals)
			colnames(vals) <- names(x)
			focals <- sp::SpatialPointsDataFrame(focals, data=vals, proj4string=crs)
		
		} else if ('SpatialPoints' %in% class(focals)) {
		
			crs <- sp::CRS(raster::projection(focals))

			vals <- raster::extract(x, focals)
			vals <- as.data.frame(vals)
			colnames(vals) <- names(x)
			focals <- sp::SpatialPointsDataFrame(focals, data=vals, crs)
			
		} else if ('SpatialPointsDataFrame' %in% class(focals)) {

			vals <- raster::extract(x, focals)
			vals <- as.data.frame(vals)
			colnames(vals) <- names(x)
			focals@data <- vals
			
		}
		
		refsAndFocalSame <- FALSE
		
	} else if (c('SpatialPointsDataFrame') %in% class(x) & is.null(focals)) {

		focals <- refs <- x
		if (is.null(vars)) vars <- names(refs)
		
		refsAndFocalSame <- TRUE
		
	} else if (any(c('matrix', 'data.frame') %in% class(x)) & is.null(focals)) {

		if (is.null(vars)) vars <- colnames(x)[3:ncol(x)]
		data <- as.data.frame(x[ , vars, drop=FALSE])
		x <- sp::SpatialPointsDataFrame(x[ , 1:2], data=data, proj4=enmSdm::getCRS('wgs84', TRUE))
		focals <- refs <- x
		
		refsAndFocalSame <- TRUE
		
	} else {
		stop('Incorrect arguments: "x" must be a raster/raster stack and "focals" NULL; OR "x" must be a raster/raster stack and "focals" a matrix, data frame, SpatialPoints, or SpatialPointsDataFrame; OR "x" must be a matrix, data frame, or SpatialPointsDataFrame, and "focals" NULL.')
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
		focals@data$DUMMY3 <- focals@data$DUMMY2 <- focals@data$DUMMY1 <- NA
		names(focals@data)[(ncol(focals@data) - 2):ncol(focals@data)] <- paste0(c('sacDistMin_', 'sacDistMid_', 'sacDistMax_'), thisVar)

		# observed values
		focalObs <- focals@data[ , thisVar, drop=TRUE]
		refsObs <- refs@data[ , thisVar, drop=TRUE]

		if (verbose) progress <- utils::txtProgressBar(min=0, max=nrow(focals), style=3, width=min(40, getOption('width')))

		# for each focals point
		for (countFocal in 1:nrow(focals)) {

			# the focals point and its data
			thisFocal <- focals[countFocal, , drop=FALSE]
						
			### distances
			
			# inter-point distances
			dists <- geosphere::distm(thisFocal, refs)
			dists <- c(dists)

			# remove distance to self
			if (refsAndFocalSame) dists[countFocal] <- NA
			
			# limit distances being considered
			if (limitDist) {
				dists <- c(dists[dists <= maxDistToInclude])
				if (length(dists) > 0) dists <- stats::na.omit(dists)
			}

			# if any inter-point distances are in the set being considered
			if (length(dists) > 0) {
				
				# distDistrib <- statisfactory::histOverlap(dists, breaks=breaks, graph=FALSE, indices=TRUE, ...)
				distDistrib <- statisfactory::histOverlap(dists, breaks=breaks, graph=FALSE, indices=TRUE)
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
				obsAbsDiffLower <- sapply(obsAbsDiff, stats::quantile, p=twoTailLower, na.rm=TRUE)
				# obsAbsDiffMid <- sapply(obsAbsDiff, quantile, twoTailMid, na.rm=TRUE)
				obsAbsDiffMid <- sapply(obsAbsDiff, mean, na.rm=TRUE)
				obsAbsDiffUpper <- sapply(obsAbsDiff, stats::quantile, p=twoTailUpper, na.rm=TRUE)

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
				
					randAbsDiffLower[iter, ] <- sapply(randAbsDiff, stats::quantile, oneTailLower, na.rm=TRUE)
					# randAbsDiffMid[iter, ] <- sapply(randAbsDiff, quantile, oneTailMid, na.rm=TRUE)
					randAbsDiffMid[iter, ] <- sapply(randAbsDiff, mean, na.rm=TRUE)
					randAbsDiffUpper[iter, ] <- sapply(randAbsDiff, stats::quantile, oneTailUpper, na.rm=TRUE)
				
				} # next iteration
			
				# calculate average lower quantile value for randomized distance distributions
				randAbsDiffLower <- apply(randAbsDiffLower, 2, mean, na.rm=TRUE)
				randAbsDiffMid <- apply(randAbsDiffMid, 2, mean, na.rm=TRUE)
				randAbsDiffUpper <- apply(randAbsDiffUpper, 2, mean, na.rm=TRUE)
		
				### remember the least/middle/most distance where the observed difference in values is no different than the randomized differences
				sacMinDistIndex <- which(randAbsDiffLower <= obsAbsDiffUpper)[1]
				
				# if there is a characteristic distance of SAC
				if (!is.na(sacMinDistIndex)) {
					focals@data[countFocal, paste0('sacDistMin_', thisVar)] <- max(0, distDistrib[sacMinDistIndex, 'lower'])
					focals@data[countFocal, paste0('sacDistMid_', thisVar)] <- max(0, distDistrib[sacMinDistIndex, 'middle'])
					focals@data[countFocal, paste0('sacDistMax_', thisVar)] <- max(0, distDistrib[sacMinDistIndex, 'upper'])
				# if none occurred in the intervals considered, return maximum value (if desired)
				} else if (returnMax) {
					numDistBins <- nrow(distDistrib)
					focals@data[countFocal, paste0('sacDistMin_', thisVar)] <- max(0, distDistrib[numDistBins, 'upper'])
					focals@data[countFocal, paste0('sacDistMid_', thisVar)] <- max(0, distDistrib[numDistBins, 'upper'])
					focals@data[countFocal, paste0('sacDistMax_', thisVar)] <- max(0, distDistrib[numDistBins, 'upper'])
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
				
			if (verbose) utils::setTxtProgressBar(progress, countFocal)
			
		} # next "focals" site
		
		if (verbose) close(progress)
		
	} # next variable
	
	# put into raster format
	if (any(c('RasterLayer', 'RasterStack', 'RasterBrick') %in% class(x)) & !focalSupplied) {
	
		for (countVar in 1:length(vars)) {
			
			focalVals <- rep(NA, raster::ncell(x))
			if (length(nonNas) > 0) {
				focalVals[nonNas] <- focals@data[ , paste0('sacDistMin_', vars[countVar])]
			} else {
				focalVals <- focals@data[ , paste0('sacDistMin_', vars[countVar])]
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
	
		focals
		
	}
	
}
