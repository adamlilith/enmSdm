#' Calculate spatial autocorrelation for geographic points
#'
#' This function calculates a measure of spatial autocorrelation for a set of geographic points. The procedure first computes the distribution of pairwise distances between the points. The frequency of distances are then tabulated in overlapping distance bins (e.g., 0 to 10 km, 5 to 15 km, 10 to 20 km, etc.). Then, many sets of randomly located points are generated. Each set contains the same number of points as observed.  For each set pairwise distances are calculated and the distance distribution tabulated. The observed distance distribution can then be compared to the randomized distributions to determine at what distance(s) the observed distances are more/less clustered than expected by chance.
#' @param x SpatialPoints or SpatialPointsDataFrame object or a matrix or data frame with two columns. If a matrix or data frame, then the first column must have values for longitude and the second latitude.
#' @param rast Raster, raster stack, or raster brick object. Randomly located sites will be placed in any non-\code{NA} cell on this raster. If this is a raster stack or brick, the first layer will be used.
#' @param breaks Three numeric values or a matrix or data frame with at least two columns:
#' \itemize{
#' \item Three numeric values: The first two values are smallest and largest distances (in km) across which to tabulate distance frequencies. The third value is the number of bins.
#' \item Matrix or data frame with at least two columns. Each row corresponds to a different bin. The first column represents the minimum distance in each bin  (in km) and the second column the maximum distance (in km). Subsequent columns are ignored. Note that by using this option arbitrary bins can be used--they need not overlap or even be continuous in coverage.
#' }
#' @param iters Positive integer, number of times to generate randomized realizations of points.
#' @param verbose Logical, if \code{TRUE} then display progress. Default is \code{FALSE}.
#' @param ... Arguments to pass to \code{\link[dismo]{randomPoints}}.
#' @details The idea behind this measure of spatial autocorrelation is that a set of geographic points is "independent" of one another if their pairwise distances are indistinguishable from pairwise distances of the same number of points randomly located across a landscape. Typically a set of points displays clustering (non-independence) across some distances but not all distances. Thus to identify the scale of clustering pairwise distances are tabulated into bins.  (We suggest using overlapping bins, e.g., from 0 to 20 km, 10 to 30 km, 20 to 40 km, etc. See the example below for how to do this).  
#' The function \code{spatialCorrForPoints} first calculates the observed distance distribution and tabulates the frequency of distances into bins. Then, it generates a set of randomly located points equal to the same number of points as in the observed set. It then calculates the randomized distance distribution and tabulates the distances. The randomization is repeated a large number of times (the default is 100). The observed frequency of distances can be compared to the set of random distances using \code{spatialCorForPointsSummary} and \code{spatialCorForPointsPlot}. The default values in those functions assume that clustering occurs if the observed pairwise distance is > the 95th quantile of the null frequency distribution for that bin (i.e., a 1-tailed test), but users can specify a different percentile to demarcate significance. In practice a series of distance bins often show clustering, but the one usually of interest is the first distance bin (the one closest to 0) that has a non-significant difference between observed and expected distances. This is the characteristic diameter of a cluster of points. Points closer than this distance can be considered non-independent of one another.
#' @seealso \code{\link[geosphere]{distm}}, \code{\link[enmSdm]{spatialCorForPointsPlot}}, \code{\link[enmSdm]{spatialCorForPointsSummary}}
#' @examples
#' \dontrun{
#' # create raster of Madagascar
#' data(mad0)
#' rast <- raster::raster(mad0, res=c(1/12, 1/12))
#' rast[] <- 1
#' rast <- raster::crop(rast, mad0)
#' mad0rast <- raster::rasterize(mad0, rast)
#' rast <- rast * mad0rast
#' 
#' # lemur point data
#' data(lemurs)
#' fulvus <- lemurs[lemurs$species == 'Eulemur fulvus', c('longitude', 'latitude')]
#' 
#' # create overlapping bins for tabulating pairwise distances
#' ext <- extent(rast)
#' southwest <- c(ext@xmin, ext@ymin)
#' northeast <- c(ext@xmax, ext@ymax)
#' maxDist <- geosphere::distGeo(southwest, northeast)
#' maxDist <- maxDist / 1000
#' 
#' binLength <- 60 # in km
#' maxDist <- binLength * ceiling(maxDist / binLength)
#' 
#' breaks <- data.frame(
#' 	lower=seq(0, maxDist - binLength, by=binLength / 2),
#' 	upper=seq(binLength, maxDist, by=binLength / 2)
#' )
#' 
#' # compare observed pairwise distance distribution to null distribution
#' # of pairwise values from randomly located points
#' obsAndNullDistrib <- spatialCorrForPoints(
#' 	x = fulvus,
#' 	rast = rast,
#' 	breaks = breaks,
#' 	iters = 100,
#' 	verbose = TRUE
#' )
#' 
#' sacDist <- spatialCorrForPointsSummary(obsAndNullDistrib)
#' main <- paste('Characteristic cluster size:', sacDist, 'km')
#' spatialCorrForPointsPlot(obsAndNullDistrib, xlab='Distance (km)', main=main)
#' }
#' @export

spatialCorrForPoints <- function(
	x,
	rast,
	breaks,
	iters = 100,
	verbose = FALSE,
	...
) {

	if (class(x) %in% 'data.frame') x <- as.matrix(x)

	if (!(class(x) %in% c('matrix', 'SpatialPoints', 'SpatialPointsDataFrame'))) {
		stop('Argument "x" must be a "SpatialPoints" or "SpatialPointsDataFrame" object or a matrix with two columns.')
	}

	numPoints <- if (class(x) == 'SpatialPoints') {
		length(x)
	} else {
		nrow(x)
	}

	# observed distances
	if (verbose) omnibus::say('Calculating observed pairwise distance distribution...')
	obsDist <- geosphere::distm(x)
	obsDist <- obsDist / 1000
	diag(obsDist) <- NA
	obsDist <- c(obsDist)
	
	distDistrib <- statisfactory::histOverlap(obsDist, breaks=breaks, graph=FALSE)
	
	# create bank of random points
	if (verbose) omnibus::say('Obtaining random points...')
	rands <- dismo::randomPoints(rast, iters * numPoints, ...)
	
	# random distributions
	if (verbose) omnibus::say('Calculating randomized distance distributions...')

	distDistrib <- distDistrib[ , c('lower', 'upper', 'proportion')]
	colnames(distDistrib)[3L] <- 'observedProportion'
	
	for (iter in 1:iters) {
	
		if (verbose) omnibus::say(iter, post=ifelse(iter %% 20 == 0 | iter == iters, 1, 0))

		randPoints <- rands[(1 + (iter - 1) * numPoints):(iter * numPoints), ]
		randDist <- geosphere::distm(randPoints)
		randDist <- randDist / 1000
		diag(randDist) <- NA
		randDist <- c(randDist)
		
		thisRandDistDistrib <- statisfactory::histOverlap(randDist, breaks=breaks, graph=FALSE)
		thisRandProportion <- thisRandDistDistrib[ , 'proportion', drop=FALSE]
		colnames(thisRandProportion) <- paste0('randProportion', iter)
		distDistrib <- cbind(distDistrib, thisRandProportion)
	
	}

	distDistrib
	
}
